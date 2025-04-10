#!/usr/bin/env python3
import pandas as pd
import sys

df = pd.read_csv(sys.argv[1], sep='\t', header=None, dtype=str)

df_unique = df.groupby(3, as_index=False).agg({0: 'first', 1: 'min', 2: 'max'})

df_unique = df_unique[[0, 1, 2, 3]]

df_unique.to_csv(sys.argv[2], sep='\t', index=False, header=False, quoting=3)

