export HF_ENDPOINT=https://hf-mirror.com
CUDA_VISIBLE_DEVICES=0 python compute_fitness.py \
    --model_name AI4Protein/ProSST-2048 \
    --model_out_name ProSST-2048 \
    --base_dir data/proteingym_v1 \
    --out_scores_dir result/proteingym_v1

export HF_ENDPOINT=https://hf-mirror.com
alpha=0.8
CUDA_VISIBLE_DEVICES=1 python compute_fitness.py \
    --model_name AI4Protein/ProSST-2048 \
    --logit_mode aa_seq_aln \
    --alpha $alpha \
    --model_out_name ProSST-2048-seq_aln"$alpha" \
    --base_dir data/proteingym_v1 \
    --out_scores_dir result/proteingym_v1

export HF_ENDPOINT=https://hf-mirror.com
alpha=0.8
CUDA_VISIBLE_DEVICES=0 python compute_fitness.py \
    --model_name AI4Protein/ProSST-2048 \
    --logit_mode struc_seq_aln \
    --alpha $alpha \
    --model_out_name ProSST-2048-struc_aln"$alpha" \
    --base_dir data/proteingym_v1 \
    --out_scores_dir result/proteingym_v1


export HF_ENDPOINT=https://hf-mirror.com
alpha=0.9
CUDA_VISIBLE_DEVICES=1 python compute_fitness.py \
    --model_name AI4Protein/ProSST-2048 \
    --logit_mode aa_seq_aln+struc_seq_aln \
    --alpha $alpha \
    --model_out_name ProSST-2048-aa_struc_aln"$alpha" \
    --base_dir data/proteingym_v1 \
    --out_scores_dir result/proteingym_v1