import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import sys
# === 修改这里：输入文件路径 ===
file_path = sys.argv[1]
output_img = 'loop_score_distribution.png'

print(f"Reading file: {file_path}")

try:
    # 读取数据 (没有表头)
    # 列定义: chr1, x1, x2, chr2, y1, y2, score
    df = pd.read_csv(file_path, sep='\t', header=None, names=['chr1', 'x1', 'x2', 'chr2', 'y1', 'y2', 'score'])
    
    # 1. 基础统计
    print("-" * 30)
    print("E Score Statistics:")
    print(df['score'].describe())
    print("-" * 30)
    
    # 2. 计算分位数 (帮助确定 cutoff)
    quantiles = df['score'].quantile([0.5, 0.75, 0.9, 0.95, 0.99])
    print("\nQuantiles (Signal Strength thresholds):")
    print(quantiles)
    
    # 3. 绘图
    plt.figure(figsize=(12, 6))
    
    # 主直方图
    sns.histplot(df['score'], bins=100, kde=True, color='skyblue', label='Distribution')
    
    # 标记常用分位线
    colors = ['green', 'orange', 'red']
    labels = ['50% (Median)', '90%', '99%']
    for q, c, l in zip([0.5, 0.9, 0.99], colors, labels):
        val = df['score'].quantile(q)
        plt.axvline(val, color=c, linestyle='--', linewidth=1.5, label=f'{l}: {val:.4f}')
    
    plt.title(f'Distribution of Loop Interaction Scores (E score)\nTotal Loops: {len(df)}')
    plt.xlabel('E Score (Interaction Strength)')
    plt.ylabel('Count')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # 保存图片
    plt.savefig(output_img, dpi=300)
    print(f"\n✅ Plot saved to: {output_img}")
    print("Tip: You can download this image to see the global distribution.")

except Exception as e:
    print(f"Error: {e}")
