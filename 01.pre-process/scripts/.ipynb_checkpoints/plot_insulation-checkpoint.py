import xarray as xr
import matplotlib.pyplot as plt
import numpy as np
import sys
import random

if len(sys.argv) < 2:
    print("Usage: python plot_insulation_random.py <path_to_nc_file>")
    sys.exit(1)

file_path = sys.argv[1]

try:
    ds = xr.open_dataset(file_path)
    n_cells = ds.dims['cell']
    cell_idx = random.randint(0, n_cells - 1)
    
    cell_name = ds.cell.values[cell_idx]
    print(f"Total cells: {n_cells}")
    print(f"Randomly selected Index: {cell_idx}")
    print(f"Cell Name: {cell_name}")

    var_name = list(ds.data_vars)[0]
    data = ds[var_name][cell_idx, :].values

    if len(data) > 600:
        plot_data = data[100:600]
        region_label = "Bins 100-600 (Approx 25Mb)"
    else:
        plot_data = data
        region_label = "All Bins"

    plt.figure(figsize=(12, 4))
    plt.plot(plot_data, label='Insulation Score', color='tab:blue', linewidth=1)
    
    plt.axhline(y=0, color='black', linewidth=0.5, alpha=0.5)

    plt.title(f"Insulation Score - {cell_name} (Index: {cell_idx})")
    plt.xlabel(f"Genomic Bins (50kb) - {region_label}")
    plt.ylabel("Score")
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # 7. 保存文件 (文件名包含细胞ID，防止覆盖)
    outfile = f"insulation_{cell_name}.png"
    plt.savefig(outfile)
    print(f"✅ Plot saved to: {outfile}")

except Exception as e:
    print(f"Error: {e}")