import os
import shutil
import argparse
from tqdm import tqdm


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='data format converter')
    parser.add_argument("--input_dir", type=str, default=None, required=True)
    parser.add_argument("--output_dir", type=str, default=None, required=True)
    parser.add_argument("--format", type=str, default="REM2SSN", choices=['SSN2REM', 'REM2SSN'], required=True)
    args = parser.parse_args()
    
    os.makedirs(args.output_dir, exist_ok=True)
    if args.format == "REM2SSN":
        output_dir = os.path.join(args.output_dir, "DATASET")
        os.makedirs(output_dir, exist_ok=True)
        proteins = os.listdir(args.input_dir + "/aa_seq")
        for p in tqdm(proteins):
            protein_name = p.split(".")[0]
            out_protein_dir = os.path.join(output_dir, protein_name)
            os.makedirs(out_protein_dir, exist_ok=True)
            shutil.copyfile(f"{args.input_dir}/aa_seq/{protein_name}.fasta", f"{out_protein_dir}/{protein_name}.fasta")
            shutil.copyfile(f"{args.input_dir}/pdbs/{protein_name}.pdb", f"{out_protein_dir}/{protein_name}.pdb")
            shutil.copyfile(f"{args.input_dir}/substitutions/{protein_name}.csv", f"{out_protein_dir}/{protein_name}.csv")
    elif args.format == "SSN2REM":
        proteins = os.listdir(args.input_dir + "/DATASET")
        pdb_dir = os.path.join(args.output_dir, "pdbs")
        os.makedirs(pdb_dir, exist_ok=True)
        aa_seq_dir = os.path.join(args.output_dir, "aa_seq")
        os.makedirs(aa_seq_dir, exist_ok=True)
        substitutions_dir = os.path.join(args.output_dir, "substitutions")
        os.makedirs(substitutions_dir, exist_ok=True)
        for p in tqdm(proteins):
            protein_name = p
            out_protein_dir = os.path.join(args.output_dir, protein_name)
            os.makedirs(out_protein_dir, exist_ok=True)
            shutil.copyfile(f"{args.input_dir}/DATASET/{protein_name}/{protein_name}.fasta", f"{out_protein_dir}/{protein_name}.fasta")
            shutil.copyfile(f"{args.input_dir}/DATASET/{protein_name}/{protein_name}.pdb", f"{out_protein_dir}/{protein_name}.pdb")
            shutil.copyfile(f"{args.input_dir}/DATASET/{protein_name}/{protein_name}.csv", f"{out_protein_dir}/{protein_name}.csv")
        