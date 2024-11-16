import argparse
import os
import pandas as pd
from tqdm import tqdm
from utils import load_coords


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='make single mutant csv')
    parser.add_argument("--fasta_dir", type=str, default=None, required=True)
    parser.add_argument("--out_dir", type=str, default=None)
    parser.add_argument("--pdb_file", type=str, default=None)
    parser.add_argument("--out_file", type=str, default=None)
    parser.add_argument("--start", type=int, default=-1)
    parser.add_argument("--end", type=int, default=int(1e6))
    args = parser.parse_args()

    one_letter = {
            'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q',
            'ASP':'D', 'ASN':'N', 'HIS':'H', 'TRP':'W', 'PHE':'F', 'TYR':'Y',
            'ARG':'R', 'LYS':'K', 'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A',
            'GLY':'G', 'PRO':'P', 'CYS':'C'
            }
    AA = list(one_letter.values())

    if args.fasta_dir is not None:
        proteins = os.listdir(args.fasta_dir)
        for p in tqdm(proteins):
            protein_name = p.split(".")[0]
            seq = open(f"{args.fasta_dir}/{p}").readlines()[1].strip()
            data = {"mutant":[], "DMS_score":[]}
            for idx, s in tqdm(enumerate(seq)):
                if idx + 1 < args.start or idx + 1 > args.end:
                    continue
                for a in AA:
                    if a == s:
                        continue
                    data["mutant"].append(f"{s}{idx+1}{a}")
                    data["DMS_score"].append(0)
            pd.DataFrame(data).to_csv(f"{args.out_dir}/{protein_name}.csv", index=False)
    
    if args.pdb_file is not None:
        out_dir = os.path.dirname(args.out_file)
        os.makedirs(out_dir, exist_ok=True)
        _, seq = load_coords(args.pdb_file, "A")
        data = {"mutant":[], "DMS_score":[]}
        for idx, s in tqdm(enumerate(seq)):
            if idx + 1 < args.start or idx + 1 > args.end:
                continue
            for a in AA:
                if a == s:
                    continue
                data["mutant"].append(f"{s}{idx+1}{a}")
                data["DMS_score"].append(0)
        pd.DataFrame(data).to_csv(args.out_file, index=False)