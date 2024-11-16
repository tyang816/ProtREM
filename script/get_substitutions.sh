protein_dir=<your_protein_dir>
python src/data/get_substitutions.py \
    --fasta_dir data/$protein_dir/aa_seq \
    --out_dir data/$protein_dir/substitutions