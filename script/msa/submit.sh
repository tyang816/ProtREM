p=BFP
d=/home/tanyang/uniref100.fasta
for b in (0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9)
do
    qsub -v protein=$p,bitscore=$b,database=$d -N "$p"_"$b" script/get_msa/get_msa.pbs
done


# evcouplings

d=/home/tanyang/data/0dataset/MEER/MEER_v1_Rep.fasta
for p in data/proteingym_v1/aa_seq/*.fasta
do
    p=$(basename $p .fasta)
    echo $p
    qsub -v protein=$p,database=$d -N $p script/get_msa/evcouplings.pbs
done


d=/home/tanyang/uniref100.fasta
for p in data/proteingym_v1/aa_seq/*.fasta
do
    p=$(basename $p .fasta)
    echo $p
    qsub -v protein=$p,database=$d -N $p script/get_msa/evcouplings.pbs
done
