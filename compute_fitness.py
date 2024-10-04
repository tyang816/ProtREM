
import torch
import os
import pandas as pd
from tqdm import tqdm
from pathlib import Path
from Bio import SeqIO
from scipy.stats import spearmanr
from transformers import AutoTokenizer, AutoModelForMaskedLM
from argparse import ArgumentParser

device = "cuda" if torch.cuda.is_available() else "cpu"


def read_multi_fasta(file_path):
    """
    params:
        file_path: path to a fasta file
    return:
        a dictionary of sequences
    """
    sequences = {}
    current_sequence = ''
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_sequence:
                    sequences[header] = current_sequence.upper()
                    current_sequence = ''
                header = line
            else:
                current_sequence += line
        if current_sequence:
            sequences[header] = current_sequence
    return sequences


def read_seq(fasta):
    for record in SeqIO.parse(fasta, "fasta"):
        return str(record.seq)

def count_matrix_from_alignment(tokenizer, alignment_file):
    alignment_dict = read_multi_fasta(alignment_file)
    aln_start, aln_end = list(alignment_dict.keys())[0].split('/')[-1].split('-')
    alignment_seqs = list(alignment_dict.values())
    print(f">>> Alignment start: {aln_start}, end: {aln_end}, start tokenizer")
    tokenized_results = tokenizer(alignment_seqs, return_tensors="pt", padding=True)
    alignment_ids = tokenized_results["input_ids"][:,1:-1]
    # count distribution of each column, [seq_len, vocab_size]
    count_matrix = torch.zeros(alignment_ids.size(1), tokenizer.vocab_size)
    for i in tqdm(range(alignment_ids.size(1))):
        count_matrix[i] = torch.bincount(alignment_ids[:,i], minlength=tokenizer.vocab_size)
    count_matrix = (count_matrix / count_matrix.sum(dim=1, keepdim=True)).to(device)
    return count_matrix, int(aln_start)-1, int(aln_end)

def tokenize_structure_sequence(structure_sequence):
    shift_structure_sequence = [i + 3 for i in structure_sequence]
    shift_structure_sequence = [1, *shift_structure_sequence, 2]
    return torch.tensor(
        [
            shift_structure_sequence,
        ],
        dtype=torch.long,
    )


@torch.no_grad()
def score_protein(model, tokenizer, residue_fasta, structure_fasta, mutant_df, alpha=0.7,
                  residu_alignment_file=None):
    sequence = read_seq(residue_fasta)
    structure_sequence = read_seq(structure_fasta)

    structure_sequence = [int(i) for i in structure_sequence.split(",")]
    ss_input_ids = tokenize_structure_sequence(structure_sequence).to(device)
    tokenized_results = tokenizer([sequence], return_tensors="pt")
    input_ids = tokenized_results["input_ids"].to(device)
    attention_mask = tokenized_results["attention_mask"].to(device)

    outputs = model(
        input_ids=input_ids,
        attention_mask=attention_mask,
        ss_input_ids=ss_input_ids,
        labels=input_ids,
    )

    logits = outputs.logits[0]
    logits = torch.log_softmax(logits[1:-1, :], dim=-1)
    
    if residu_alignment_file is not None:
        print(">>> Using residue sequence alignment matrix...")
        alignment_matrix, aln_start, aln_end = count_matrix_from_alignment(tokenizer, residu_alignment_file)
        aln_modify_logits = (1-alpha) * logits[aln_start: aln_end, :] + alpha * torch.log_softmax(alignment_matrix, dim=-1)
        logits = torch.cat([logits[:aln_start], aln_modify_logits, logits[aln_end:]], dim=0)

    mutants = mutant_df["mutant"].tolist()
    scores = []
    vocab = tokenizer.get_vocab()
    print(">>> Scoring mutants...")
    for mutant in tqdm(mutants):
        pred_score = 0
        for sub_mutant in mutant.split(":"):
            wt, idx, mt = sub_mutant[0], int(sub_mutant[1:-1]) - 1, sub_mutant[-1]
            score = logits[idx, vocab[mt]] - logits[idx, vocab[wt]]
            pred_score += score.item()
        scores.append(pred_score)

    return scores
    


def read_names(fasta_dir):
    files = Path(fasta_dir).glob("*.fasta")
    names = [file.stem for file in files]
    return names


def main():
    parser = ArgumentParser()
    parser.add_argument("--model_name", type=str, default="AI4Protein/ProSST-2048", nargs="+", required=True)
    parser.add_argument("--model_out_name", type=str, default="ProSST-2048", nargs="+", help="Output model name",)
    
    # data directories
    parser.add_argument("--base_dir", type=str, default=None, help="Base directory containing all data",)
    parser.add_argument("--residue_dir", type=str, default=None, help="Directory containing FASTA files of residue sequences",)
    parser.add_argument("--structure_dir", type=str, default=None, help="Directory containing FASTA files of structure sequences",)
    parser.add_argument("--mutant_dir", type=str, default=None, help="Directory containing CSV files with mutants",)
    
    # retrieval and logits mode
    parser.add_argument("--retrieval_mode", type=str, default=None, choices=["seq", "struc", "seq+struc"], help="Mode to retrieve data",)
    parser.add_argument("--alpha", type=float, default=0.5, help="Alpha value for combining logits",)
    parser.add_argument("--residue_alignment_dir", type=str, default=None, help="Directory containing a2m files of residue alignments",)
    
    # output directory
    parser.add_argument("--out_scores_dir", default=None, help="Directory to save scores")
    args = parser.parse_args()

    print("Scoring proteins...")
    os.makedirs(args.out_scores_dir, exist_ok=True)
    os.makedirs(f"{args.out_scores_dir}/scores", exist_ok=True)
    if args.base_dir:
        protein_names = read_names(f"{args.base_dir}/residue_sequence")
    if args.residue_dir:
        protein_names = read_names(args.residue_dir)
    protein_names = sorted(protein_names)
    corrs = []
    print(protein_names)
    print(f">>> total proteins: {len(protein_names)}")
    print("=====================================")
    for model_idx, model_name in enumerate(args.model_name):
        print(f">>> Loading model {model_name}...")
        model = AutoModelForMaskedLM.from_pretrained(
            model_name, trust_remote_code=True
        )
        model = model.to(device)
        tokenizer = AutoTokenizer.from_pretrained(model_name, trust_remote_code=True)
        
        for idx, protein_name in enumerate(protein_names):
            print(f">>> Scoring {protein_name}, current {idx+1}/{len(protein_names)}...")
            # load data
            if args.base_dir:
                residue_fasta = f"{args.base_dir}/residue_sequence/{protein_name}.fasta"
                structure_fasta = f"{args.base_dir}/structure_sequence/{model_name.split('-')[-1]}/{protein_name}.fasta"
                mutant_file = f"{args.base_dir}/substitutions/{protein_name}.csv"
                if args.retrieval_mode == "seq":
                    residu_alignment_file = f"{args.base_dir}/residue_alignment/{protein_name}.a2m"
                else:
                    residu_alignment_file = None
            if args.residue_dir:
                residue_fasta = f"{args.residue_dir}/{protein_name}.fasta"
            if args.structure_dir:
                structure_fasta = f"{args.structure_dir}/{protein_name}.fasta"
            if args.mutant_dir:
                mutant_file = f"{args.mutant_dir}/{protein_name}.csv"
            if os.path.exists(f"{args.out_scores_dir}/scores/{protein_name}.csv"):
                mutant_file = f"{args.out_scores_dir}/scores/{protein_name}.csv"
            mutant_df = pd.read_csv(mutant_file)
            
            if args.model_out_name:
                model_out_name = args.model_out_name[model_idx]
            else:
                model_out_name = model_name.split("/")[-1]
                
            if model_out_name not in mutant_df.columns:
                scores = score_protein(
                        model=model,
                        tokenizer=tokenizer,
                        residue_fasta=residue_fasta,
                        structure_fasta=structure_fasta,
                        mutant_df=mutant_df,
                        alpha=args.alpha,
                        residu_alignment_file=residu_alignment_file,
                    )
                
                mutant_df[model_out_name] = scores
            
            corr = spearmanr(mutant_df["DMS_score"], mutant_df[model_out_name]).correlation
            corrs.append(corr)
            print(f">>> {model_name} on {protein_name}: {corr}")
            print("=====================================")
            mutant_df.to_csv(f"{args.out_scores_dir}/scores/{protein_name}.csv", index=False)
        print(f"====== {model_name} average correlation: {sum(corrs)/len(corrs)} ======")
        pd.DataFrame({'protein': protein_names, model_out_name: corrs}).to_csv(f"{args.out_scores_dir}/correlation.csv", index=False)

if __name__ == "__main__":
    main()