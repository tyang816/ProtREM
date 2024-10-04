export HF_ENDPOINT=https://hf-mirror.com
CUDA_VISIBLE_DEVICES=0 python compute_fitness.py \
    --model_name AI4Protein/ProSST-2048 \
    --model_out_name ProSST-2048 \
    --base_dir data/proteingym_v1 \
    --out_scores_dir result/proteingym_v1


export HF_ENDPOINT=https://hf-mirror.com
CUDA_VISIBLE_DEVICES=3 python compute_fitness.py \
    --model_name AI4Protein/ProSST-2048 \
    --retrieval_mode seq \
    --alpha 0.3 \
    --model_out_name ProSST-2048-seq_aln0.3 \
    --base_dir data/proteingym_v1 \
    --out_scores_dir result/proteingym_v1