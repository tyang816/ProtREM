raw_dir=data/proteingym_v1
# loop the files of data/proteingym_v1/aa_seq_aln_a3m_raw
for file in $raw_dir/aa_seq_aln_a3m_raw/*
do
    file=$(basename $file)
    reformat.pl $raw_dir/aa_seq_aln_a3m_raw/$file $raw_dir/aa_seq_aln_a3m/$file -r
done