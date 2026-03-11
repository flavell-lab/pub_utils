"""
Compute network metrics for all preassembled connectomes.

Outputs: connectomes/metrics/preassembled_metrics.csv

Usage:
    python scripts/compute_preassembled_metrics.py
"""

from pub_utils.metrics import compute_all_metrics

summary = compute_all_metrics(save_path="connectomes/metrics/preassembled_metrics.csv")
print(f"\nDone. {len(summary)} connectomes processed.")
print(summary[["label", "density", "Q", "reciprocity", "total_degree_cv"]].to_string(index=False))
