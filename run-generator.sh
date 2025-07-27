#!/bin/bash

# === 参数说明 ===
# $1 = dataset（如 enron-email）
# $2 = k
# $3 = t
# $4 = max
# $5 = min

# === usage 信息 ===
echo_usage() {
    echo
    echo "usage: $(basename "$0") [dataset] [k-limit] [thread] [max-weight] [min-weight]"
    echo "e.g  : $(basename "$0") enron-email 3 8 100 1"
    echo
    echo "dataset supported:"
    echo "   - enron-email"
    echo "   - twitch"
    echo "   - ...."
    echo
}

BIN_DIR="$(dirname "$0")/bin"
DATA_DIR="$(dirname "$0")/data"
EXEC="$BIN_DIR/experiment_program"

# 检查参数
if [ -z "$1" ]; then
    echo "[ERROR] Missing parameter: dataset"
    echo_usage
    exit 1
fi

if [ -z "$2" ]; then
    echo "[ERROR] Missing parameter: k"
    echo_usage
    exit 1
fi

if [ -z "$3" ]; then
    echo "[ERROR] Missing parameter: t"
    echo_usage
    exit 1
fi

if [ -z "$4" ]; then
    echo "[ERROR] Missing parameter: max"
    echo_usage
    exit 1
fi

if [ -z "$5" ]; then
    echo "[ERROR] Missing parameter: min"
    echo_usage
    exit 1
fi

dataset="$1"
k="$2"
t="$3"
max="$4"
min="$5"

# 根据 dataset 设置路径
case "$dataset" in
    enron-email)
        file="$DATA_DIR/enron-email/raw/Email-Enron.txt"
        out="$DATA_DIR/enron-email/processed/"
        ;;
    twitch)
        file="$DATA_DIR/twitch/raw/large_twitch_edges.txt"
        out="$DATA_DIR/twitch/processed/"
        ;;
    *)
        echo "[ERROR] unknown dataset: $dataset"
        echo_usage
        exit 1
        ;;
esac

echo
echo "[INFO] Running label generation..."
echo "$EXEC generate-label -t $t -f $file -p $out -k $k -max $max -min $min"
"$EXEC" generate-label -t "$t" -f "$file" -p "$out" -k "$k" -max "$max" -min "$min"

echo
echo "[INFO] Plotting degree distribution with Python..."
python3 "$BIN_DIR/plot_degree.py" "$out"

exit 0


