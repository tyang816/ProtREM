# 查询序列文件
QUERY_FILE="data/case/aa_seq/BFP.fasta"
# 目标数据库
DATABASE="/home/tanyang/data/0dataset/MEER/MEER_0.5.fasta"
# 输出文件夹
OUTPUT_DIR="jackhmmer_results"

# 创建输出文件夹
mkdir -p $OUTPUT_DIR

# 设置迭代次数
ITERATIONS=5

# 比特得分阈值范围（从0.1到0.9，共9个阈值）
BIT_SCORES=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9)

for BIT_SCORE in "${BIT_SCORES[@]}"; do
    # 输出文件前缀
    OUTPUT_PREFIX="${OUTPUT_DIR}/results_bitscore_${BIT_SCORE}"

    # 运行Jackhmmer，保存对齐结果为Stockholm格式
    jackhmmer -N "$ITERATIONS" --incT "$BIT_SCORE" -A "${OUTPUT_PREFIX}.sto" \
        --tblout "${OUTPUT_PREFIX}_table.txt" "$QUERY_FILE" "$DATABASE" > "${OUTPUT_PREFIX}_output.txt"

    # 将Stockholm格式转换为A2M格式
    esl-reformat a2m "${OUTPUT_PREFIX}.sto" > "${OUTPUT_PREFIX}.a2m"

    echo "完成比特得分阈值为 ${BIT_SCORE} 的Jackhmmer搜索，并保存为A2M格式"
done

# 初始化变量，用于选择最佳MSA
BEST_MSA=""
MAX_EC_COUNT=0

# 遍历每个生成的MSA，计算显著的ECs数量
for BIT_SCORE in "${BIT_SCORES[@]}"; do
    OUTPUT_PREFIX="${OUTPUT_DIR}/results_bitscore_${BIT_SCORE}"
    ALIGNMENT_FILE="${OUTPUT_PREFIX}.a2m"

    # 使用ECs计算工具（例如EVcouplings）计算进化耦合
    # 假设计算结果保存在 "${OUTPUT_PREFIX}_ecs.txt"
    # 这里是占位符命令，需要替换为实际的EC计算命令
    evcouplings_cli protein --msa "${ALIGNMENT_FILE}" --other_parameters > "${OUTPUT_PREFIX}_ecs.txt"

    # 统计显著的ECs数量（根据实际的显著性阈值和输出格式调整命令）
    SIGNIFICANT_EC_COUNT=$(grep -c "significant" "${OUTPUT_PREFIX}_ecs.txt")

    echo "比特得分 ${BIT_SCORE}: ${SIGNIFICANT_EC_COUNT} 个显著的ECs"

    # 更新最佳MSA
    if [ "$SIGNIFICANT_EC_COUNT" -gt "$MAX_EC_COUNT" ]; then
        MAX_EC_COUNT="$SIGNIFICANT_EC_COUNT"
        BEST_MSA="$ALIGNMENT_FILE"
    fi
done

echo "最佳MSA是 ${BEST_MSA}，具有 ${MAX_EC_COUNT} 个显著的ECs"