protein_dir=<your_protein_dir>
output_dir=<your_output_dir>
python src/data/data_format_convert.py \
    --input_dir data/$protein_dir \
    --output_dir $output_dir \
    --format REM2SSN