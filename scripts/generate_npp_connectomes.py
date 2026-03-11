"""
Generate all 92 individual NPP-GPCR mid-range connectomes as standardized 302x302 CSVs.

Reads raw data from data/RipollSanchez2023/individual/individual_mid_range/{001..092}.csv
and writes standardized CSVs to connectomes/preassembled/molecular/{NPP}_{GPCR}_RipollSanchez2023_midRange.csv

Usage:
    python scripts/generate_npp_connectomes.py
"""

import pandas as pd
from pub_utils.io import create_neuropeptide_mapping, standardize_dataframe
from pub_utils.constants import AllHermNeurons
from pathlib import Path

RAW_DIR = Path("data/RipollSanchez2023/individual/individual_mid_range")
OUT_DIR = Path("connectomes/preassembled/molecular")

mapping = create_neuropeptide_mapping()

generated = 0
skipped = 0

for pair_key, file_idx in sorted(mapping.items(), key=lambda x: x[1]):
    # pair_key is lowercase "flp-10 dmsr-1"; convert to uppercase for filename
    npp, gpcr = pair_key.upper().split(" ")
    out_name = f"{npp}_{gpcr}_RipollSanchez2023_midRange.csv"
    out_path = OUT_DIR / out_name

    if out_path.exists():
        skipped += 1
        continue

    raw_path = RAW_DIR / f"{file_idx}.csv"
    df = pd.read_csv(raw_path)
    std_df = standardize_dataframe(df, AllHermNeurons, verbose=False)
    std_df.to_csv(out_path)
    generated += 1
    print(f"  [{file_idx}] {out_name}")

print(f"\nDone. Generated {generated}, skipped {skipped} (already existed).")
