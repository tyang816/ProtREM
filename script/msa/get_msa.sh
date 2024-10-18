database=~/uniref50.fasta
protein=BFP
bitscore=(0.2 0.3)
python src/data/get_msa.py \
    --query_file data/case/aa_seq/"$protein".fasta \
    --bitscores $bitscore \
    --database $database \
    --output_dir jackhmmer_results/$protein/"$protein"_uniref100