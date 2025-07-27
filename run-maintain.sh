#!/bin/bash

# 参数说明
# $1 = dataset（如 enron-email）
# $2 = k
# $3 = t
# $4 = max
# $5 = min
# $6 = m 迭代次数
# $7 = c 变化次数
# $8 = strategy 变化策略
# $9 = (可选) check

# ==== usage 函数（必须提前定义）====
echo_usage() {
    echo
    echo "usage: $(basename "$0") [dataset] [k-limit] [thread] [max-weight] [min-weight] [iteration_count] [edge_change_count] [strategy] [check]"
    echo "e.g  : $(basename "$0") enron-email 3 8 100 1 10 500 high_high_increase check"
    echo
    echo "dataset supported:"
    echo "   - enron-email"
    echo "strategy supported:"
    echo "   - high_high_increase"
    echo "   - high_high_decrease"
    echo "   - high_high_mixed"
    echo "   - high_low_increase"
    echo "   - high_low_decrease"
    echo "   - high_low_mixed"
    echo "   - low_low_increase"
    echo "   - low_low_decrease"
    echo "   - low_low_mixed"
}

# ==== 参数校验 ====
if [ -z "$1" ]; then echo "[ERROR] Missing parameter: dataset"; echo_usage; exit 1; fi
if [ -z "$2" ]; then echo "[ERROR] Missing parameter: k"; echo_usage; exit 1; fi
if [ -z "$3" ]; then echo "[ERROR] Missing parameter: t"; echo_usage; exit 1; fi
if [ -z "$4" ]; then echo "[ERROR] Missing parameter: max"; echo_usage; exit 1; fi
if [ -z "$5" ]; then echo "[ERROR] Missing parameter: min"; echo_usage; exit 1; fi
if [ -z "$6" ]; then echo "[ERROR] Missing parameter: m"; echo_usage; exit 1; fi
if [ -z "$7" ]; then echo "[ERROR] Missing parameter: c"; echo_usage; exit 1; fi
if [ -z "$8" ]; then echo "[ERROR] Missing parameter: strategy"; echo_usage; exit 1; fi

# ==== 变量赋值 ====
dataset="$1"
k="$2"
t="$3"
max="$4"
min="$5"
m="$6"
c="$7"
strategy="$8"
check="$9"

BIN_DIR="$(dirname "$0")/bin"
DATA_DIR="$(dirname "$0")/data"
EXEC="$BIN_DIR/experiment_program"  # Linux 下建议使用无 .exe 后缀

# ==== 支持的策略列表 ====
STRATEGIES=(
    high_high_increase high_high_decrease high_high_mixed
    high_low_increase high_low_decrease high_low_mixed
    low_low_increase low_low_decrease low_low_mixed
)

# ==== 检查 strategy 合法性 ====
found=0
for s in "${STRATEGIES[@]}"; do
    if [[ "$strategy" == "$s" ]]; then
        found=1
        break
    fi
done

if [ "$found" -ne 1 ]; then
    echo "[ERROR] Invalid strategy: $strategy"
    echo -n "Valid options: "
    printf "%s " "${STRATEGIES[@]}"
    echo
    echo_usage
    exit 1
fi

# ==== 可选参数 check ====
if [[ "$check" == "check" ]]; then
    check="--$check"
else
    check=""
fi

# ==== 数据集路径 ====
case "$dataset" in
    enron-email)
        file="$DATA_DIR/enron-email/processed/"
        out="$DATA_DIR/enron-email/processed/"
        ;;
    twitch)
        file="$DATA_DIR/twitch/processed/"
        out="$DATA_DIR/twitch/processed/"
        ;;
    *)
        echo "[ERROR] Unknown dataset: $dataset"
        echo_usage
        exit 1
        ;;
esac

# ==== 执行命令 ====
echo
echo "[INFO] Running label maintain..."
echo "$EXEC maintain-label -t $t -f $file -p $out -k $k -c $c -m $m -max $max -min $min --strategy $strategy $check -n $dataset"
for ((i = 1; i <= 10; i++)); do
"$EXEC" maintain-label -t "$t" -f "$file" -p "$out" -k "$k" -c "$c" -m "$m" -max "$max" -min "$min" --strategy "$strategy" $check -n "$dataset"
done
