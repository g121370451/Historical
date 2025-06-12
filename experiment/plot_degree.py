# import pandas as pd
# import matplotlib.pyplot as plt
# import numpy as np
# import sys
# import os

# def main():
#     if len(sys.argv) < 2:
#         print("Usage: python plot_degree.py <output_directory>")
#         sys.exit(1)

#     output_dir = sys.argv[1]
#     # output_dir = "D:\\project\\mine\\test_cmake\\experiment\\email-enron"
#     os.makedirs(output_dir, exist_ok=True)

#     input_path = os.path.join(output_dir, "degree_dist.csv")
#     output_path = os.path.join(output_dir, "degree_distribution.png")

#     df = pd.read_csv(input_path)

#     degrees = df["Degree"]
#     degrees_nonzero = degrees[degrees > 0]  # 过滤掉零度节点

#     log_degrees = np.log2(degrees_nonzero)
#     log_low = np.percentile(log_degrees, 33)
#     log_high = np.percentile(log_degrees, 67)

#     # 转回原始度数
#     low_thresh = int(np.floor(2 ** log_low))
#     high_thresh = int(np.ceil(2 ** log_high))

#     print(f"Low threshold: <= {low_thresh}, High threshold: >= {high_thresh}")

#     # 分段数据
#     low_part = df[df["Degree"] <= low_thresh]
#     mid_part = df[(df["Degree"] > low_thresh) & (df["Degree"] < high_thresh)]
#     high_part = df[df["Degree"] >= high_thresh]

#     # 画图
#     plt.figure(figsize=(10, 6))

#     # 主曲线分段（三段不同颜色）
#     plt.plot(low_part["Degree"], low_part["Count"], label=f"Low (≤{low_thresh})", color="green")
#     plt.plot(mid_part["Degree"], mid_part["Count"], label="Mid", color="blue")
#     plt.plot(high_part["Degree"], high_part["Count"], label=f"High (≥{high_thresh})", color="red")

#     # 设置图形属性
#     plt.xlabel("Degree")
#     plt.ylabel("Count")
#     plt.title("Degree Distribution")
#     plt.yscale("log")
#     # plt.xscale("log")  # 可选：若要扩展高维度可启用
#     plt.grid(True, linestyle="--", alpha=0.6)
#     plt.legend()
#     plt.tight_layout()

#     # 保存图像
#     plt.savefig(output_path)
#     print(f"Plot saved to {output_path}")

# if __name__ == "__main__":
#     main()
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os

def load_and_group_by_degree(csv_path):
    df = pd.read_csv(csv_path)
    grouped = df.groupby("Degree").size().reset_index(name="Count")
    return df, grouped

def main():
    if len(sys.argv) < 2:
        print("Usage: python plot_degree.py <output_directory>")
        sys.exit(1)

    output_dir = sys.argv[1]
    # output_dir = "D:\\project\\mine\\test_cmake\\experiment\\email-enron"
    os.makedirs(output_dir, exist_ok=True)
    # 假设这三个 CSV 在当前目录
    #     input_path = os.path.join(output_dir, "degree_dist.csv")
    #     output_path = os.path.join(output_dir, "degree_distribution.png")
    low_path = os.path.join(output_dir,  "low.csv")
    mid_path = os.path.join(output_dir, "mid.csv")
    high_path = os.path.join(output_dir,  "high.csv")
    output_path = os.path.join(output_dir,  "degree_distribution.png")

    # 读取 & groupBy degree
    low_df, low_group = load_and_group_by_degree(low_path)
    mid_df, mid_group = load_and_group_by_degree(mid_path)
    high_df, high_group = load_and_group_by_degree(high_path)

    # 找最大 vertex ID 所在的度数
    split1 = mid_df.loc[mid_df["Vertex"].idxmax()]["Degree"]
    split1Vertex = mid_df.loc[mid_df["Vertex"].idxmax()]["Vertex"]
    split2 = high_df.loc[high_df["Vertex"].idxmax()]["Degree"]
    split2Vertex = high_df.loc[high_df["Vertex"].idxmax()]["Vertex"]
    split3 = low_df.loc[low_df["Vertex"].idxmax()]["Degree"]
    split3Vertex = low_df.loc[low_df["Vertex"].idxmax()]["Vertex"]

    # 绘图
    plt.figure(figsize=(10, 6))
    plt.plot(low_group["Degree"], low_group["Count"], label="Low Degree", color="green", marker="o")
    plt.plot(mid_group["Degree"], mid_group["Count"], label="Mid Degree", color="blue", marker="^")
    plt.plot(high_group["Degree"], high_group["Count"], label="High Degree", color="orange", marker="s")

    # 标记两条竖线
    plt.axvline(x=split1, color='purple', linestyle='--', label=f"Split Mid: {split1Vertex}")
    plt.axvline(x=split2, color='red', linestyle='--', label=f"Split High: {split2Vertex}")
    plt.axvline(x=split3, color='pink', linestyle='--', label=f"Split High: {split3Vertex}")

    # 图形美化
    plt.xlabel("Degree")
    plt.ylabel("Count")
    plt.title("Degree Distribution by Class")
    plt.yscale("log")
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend()
    plt.tight_layout()

    # 保存
    plt.savefig(output_path)
    print("Saved to degree_distribution_split.png")

if __name__ == "__main__":
    main()
