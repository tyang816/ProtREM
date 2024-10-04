import os
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_dir", type=str, default=None, help="Directory containing PDB files",)
    parser.add_argument("--out_dir", type=str, default=None, help="Output directory",)