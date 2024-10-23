protein=OFP
database=/home/tanyang/uniref100.fasta
evcouplings \
    -P output/case/$protein/$protein \
    -p $protein \
    -s data/case/aa_seq/$protein.fasta \
    -d $database \
    -b "0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9" \
    -n 5 src/single_config_monomer.txt


protein=VHH_tm
database=/home/tanyang/data/0dataset/MEER/MEER_v1_Rep.fasta
evcouplings \
    -P output_MEER/case/$protein/$protein \
    -p $protein \
    -s data/case/aa_seq/$protein.fasta \
    -d $database \
    -b "0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9" \
    -n 5 src/single_config_monomer.txt