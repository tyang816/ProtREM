protein=BFP
database=/home/tanyang/uniref50.fasta
evcouplings \
    -P output/$protein/$protein \
    -p $protein \
    -s /home/tanyang/workspace/ProSST-R/data/case/aa_seq/$protein.fasta \
    -d $database \
    -b "0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9" \
    -n 5 src/single_config_monomer.txt