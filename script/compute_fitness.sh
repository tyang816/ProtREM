export HF_ENDPOINT=https://hf-mirror.com
CUDA_VISIBLE_DEVICES=0 python compute_fitness.py \
    --model_name AI4Protein/ProSST-2048 \
    --model_out_name ProSST-2048 \
    --base_dir data/proteingym_v1 \
    --out_scores_dir result/proteingym_v1

export HF_ENDPOINT=https://hf-mirror.com
CUDA_VISIBLE_DEVICES=0 python compute_fitness.py \
    --model_name AI4Protein/ProSST-2048 \
    --logit_mode aa_seq_aln \
    --alpha 0.9 \
    --model_out_name ProSST-2048-seq_aln0.9 \
    --base_dir data/proteingym_v1 \
    --out_scores_dir result/proteingym_v1

export HF_ENDPOINT=https://hf-mirror.com
CUDA_VISIBLE_DEVICES=1 python compute_fitness.py \
    --model_name AI4Protein/ProSST-2048 \
    --logit_mode struc_seq_aln \
    --alpha 0.5 \
    --model_out_name ProSST-2048-struc_aln0.5 \
    --base_dir data/proteingym_v1 \
    --out_scores_dir result/proteingym_v1