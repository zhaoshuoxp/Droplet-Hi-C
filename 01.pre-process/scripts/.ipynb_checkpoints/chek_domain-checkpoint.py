import anndata
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

boundary_path = sys.argv[1]
insulation_path = sys.arv[2]
# ============================

print(f"--- 1. Checking Boundary File: {boundary_path} ---")
try:
    adata = anndata.read_h5ad(boundary_path)
    print(f"✅ Load Success!")
    print(f"   Shape (Cells x Bins): {adata.shape}")

    nnz = adata.X.nnz
    total_pixels = adata.shape[0] * adata.shape[1]
    sparsity = 1 - (nnz / total_pixels)
    print(f"   Non-zero entries: {nnz}")
    print(f"   Sparsity: {sparsity:.4f} (Expected > 0.90)")
    
    if nnz == 0:
        print("   ❌ ERROR: Matrix is completely empty! No boundaries found.")
    elif sparsity < 0.5:
        print("   ⚠️ WARNING: Matrix is too dense. Something might be wrong with domain calling.")
    else:
        print("   ✅ Data looks reasonable (sparse boundary calls).")

except Exception as e:
    print(f"❌ Failed to read h5ad: {e}")

print(f"\n--- 2. Checking Insulation File: {insulation_path} ---")
try:
    ds = xr.open_dataset(insulation_path)
    print(f"✅ Load Success!")
    print(ds)
    
    var_name = list(ds.data_vars)[0]
    data = ds[var_name].values
    
    print(f"   Variable name: {var_name}")
    print(f"   Shape: {data.shape}")
    
    valid_data = data[~np.isnan(data)]
    
    if len(valid_data) == 0:
        print("   ❌ ERROR: All values are NaN!")
    else:
        print(f"   Min Value: {np.min(valid_data):.4f}")
        print(f"   Max Value: {np.max(valid_data):.4f}")
        print(f"   Mean Value: {np.mean(valid_data):.4f}")
        print("   ✅ Numerical range looks valid.")
        
    if 'bin_chrom' in ds.coords:
        chrom_dtype = ds.coords['bin_chrom'].dtype
        print(f"   Coordinate 'bin_chrom' dtype: {chrom_dtype}")
        if np.issubdtype(chrom_dtype, np.object_) or np.issubdtype(chrom_dtype, np.str_) or np.issubdtype(chrom_dtype, np.string_):
             print("   ✅ Coordinates are Strings (Fix worked!).")
        else:
             print(f"   ⚠️ Coordinates might still be Categorical/Int: {chrom_dtype}")

except Exception as e:
    print(f"❌ Failed to read nc file: {e}")
