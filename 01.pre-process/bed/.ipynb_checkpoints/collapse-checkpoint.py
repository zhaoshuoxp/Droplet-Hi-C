#!/usr/bin/env python3
import pandas as pd
import sys

# 读取 TSV 文件，确保所有列为字符串格式
df = pd.read_csv(sys.argv[1], sep='\t', header=None, dtype=str)

# 对第四列去重，同时第二列取最小值，第三列取最大值
df_unique = df.groupby(3, as_index=False).agg({0: 'first', 1: 'min', 2: 'max'})

# 保持原始列顺序
df_unique = df_unique[[0, 1, 2, 3]]

# 保存为 TSV，保持原始格式
df_unique.to_csv(sys.argv[2], sep='\t', index=False, header=False, quoting=3)

