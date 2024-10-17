p=BFP
d=/home/tanyang/uniref100.fasta
for b in (0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9)
do
    qsub -v protein=$p,bitscore=$b,database=$d -N "$p"_"$b" script/get_msa/get_msa.pbs
done