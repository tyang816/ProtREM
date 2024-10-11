# defalut setting
export HF_ENDPOINT=https://hf-mirror.com
CUDA_VISIBLE_DEVICES=0 python compute_fitness.py \
    --model_name AI4Protein/ProSST-2048 \
    --model_out_name ProSST-2048 \
    --base_dir data/proteingym_v1 \
    --out_scores_dir result/proteingym_v1


# zero-shot with aa sequence alignment
export HF_ENDPOINT=https://hf-mirror.com
for alpha in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
do
    CUDA_VISIBLE_DEVICES=0 python compute_fitness.py \
        --model_name AI4Protein/ProSST-2048 \
        --logit_mode aa_seq_aln \
        --alpha $alpha \
        --model_out_name ProSST-2048-seq_aln"$alpha" \
        --base_dir data/proteingym_v1 \
        --out_scores_dir result/proteingym_v1
done

export HF_ENDPOINT=https://hf-mirror.com
alpha=0.8
CUDA_VISIBLE_DEVICES=0 python compute_fitness.py \
    --model_name AI4Protein/ProSST-2048 \
    --logit_mode aa_seq_aln \
    --alpha $alpha \
    --model_out_name ProSST-2048-seq_aln"$alpha"-a3m \
    --base_dir data/proteingym_v1 \
    --aa_seq_aln_dir aa_seq_aln_a3m \
    --out_scores_dir result/proteingym_v1


# zero-shot with structure sequence alignment
export HF_ENDPOINT=https://hf-mirror.com
for alpha in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
do
    CUDA_VISIBLE_DEVICES=0 python compute_fitness.py \
        --model_name AI4Protein/ProSST-2048 \
        --logit_mode struc_seq_aln \
        --alpha $alpha \
        --model_out_name ProSST-2048-struc_aln"$alpha" \
        --base_dir data/proteingym_v1 \
        --out_scores_dir result/proteingym_v1
done

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

export HF_ENDPOINT=https://hf-mirror.com
alpha=0.8
CUDA_VISIBLE_DEVICES=1 python compute_fitness.py \
    --model_name AI4Protein/ProSST-2048 \
    --logit_mode struc_seq_aln+aa_seq_aln \
    --alpha $alpha \
    --model_out_name ProSST-2048-aa_struc_aln"$alpha" \
    --base_dir data/proteingym_v1 \
    --out_scores_dir result/proteingym_v1