"""
Visualize network metrics for preassembled connectomes.

Reads:  connectomes/metrics/preassembled_metrics.csv
Saves:  connectomes/metrics/density.png
        connectomes/metrics/modularity.png
        connectomes/metrics/reciprocity.png
        connectomes/metrics/hubness.png
        connectomes/metrics/clustering.png
        connectomes/metrics/density_clean.png
        connectomes/metrics/modularity_clean.png
        connectomes/metrics/reciprocity_clean.png
        connectomes/metrics/hubness_clean.png
        connectomes/metrics/clustering_clean.png
        connectomes/metrics/preassembled_metrics_heatmap.png

Usage:
    python scripts/plot_metrics.py
"""

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import numpy as np
from adjustText import adjust_text

# ── Load data ────────────────────────────────────────────────────────────────
df = pd.read_csv("connectomes/metrics/preassembled_metrics.csv")

# Short display name: just the dataset column
df["name"] = df["dataset"]

# Consistent type ordering and colors
TYPE_ORDER = ["chemical_synapse", "electrical_synapse", "neuropeptide", "monoamine"]
TYPE_LABELS = ["Chemical\nsynapse", "Electrical\nsynapse", "Neuropeptide", "Monoamine"]
TYPE_COLORS = ["#E24A33", "#348ABD", "#988ED5", "#8EBA42"]
type_color = dict(zip(TYPE_ORDER, TYPE_COLORS))
type_label = dict(zip(TYPE_ORDER, TYPE_LABELS))

df["color"] = df["connectome_type"].map(type_color)
df["type_label"] = df["connectome_type"].map(type_label)

# Sort so types are grouped
df["type_rank"] = df["connectome_type"].map({t: i for i, t in enumerate(TYPE_ORDER)})
df = df.sort_values(["type_rank", "density"], ascending=[True, False]).reset_index(drop=True)

OUT_DIR = "connectomes/metrics"


# ── Top/bottom annotation per group ─────────────────────────────────────────
def get_extreme_mask(series, top_n=2, bot_n=2):
    """Return boolean Series selecting the top and bottom values."""
    if len(series) <= top_n + bot_n:
        return pd.Series(True, index=series.index)
    idx = series.nlargest(top_n).index
    if bot_n > 0:
        idx = idx.union(series.nsmallest(bot_n).index)
    return pd.Series(series.index.isin(idx), index=series.index)


# ── Individual strip-plot figures ────────────────────────────────────────────
METRICS = [
    ("density",                "Edge Density",                True,  "density"),
    ("Q",                      "Modularity (Q)",              False, "modularity"),
    ("reciprocity",            "Reciprocity",                 False, "reciprocity"),
    ("total_degree_cv",        "Total Degree CV (hubness)",   True,  "hubness"),
    ("clustering_coefficient", "Clustering Coefficient",      False, "clustering"),
]


def make_strip_plot(col, ylabel, log_scale, filename, annotate=True, top_n=2, bot_n=2, ylim_bottom=None):
    """Draw a single strip-plot metric figure and save it."""
    fig, ax = plt.subplots(figsize=(6, 6))
    rng = np.random.default_rng(42)

    texts = []  # for adjust_text

    for ctype in TYPE_ORDER:
        sub = df[df["connectome_type"] == ctype].copy()
        jitter = rng.uniform(-0.08, 0.08, len(sub))
        x_pos = TYPE_ORDER.index(ctype)

        ax.scatter(
            x_pos + jitter,
            sub[col],
            c=type_color[ctype],
            s=40,
            alpha=0.75,
            edgecolors="white",
            linewidths=0.4,
            zorder=3,
        )

        # Median bar
        median = sub[col].median()
        ax.hlines(median, x_pos - 0.3, x_pos + 0.3, color=type_color[ctype],
                  linewidth=2.5, zorder=4)

        # Annotate top-2 / bottom-2
        if annotate:
            extreme_mask = get_extreme_mask(sub[col], top_n=top_n, bot_n=bot_n)
            for idx, (_, row) in enumerate(sub.iterrows()):
                if extreme_mask.iloc[idx]:
                    label_text = row['name']
                    if ctype == "neuropeptide" and "_RipollSanchez2023_midRange" in label_text:
                        label_text = label_text.replace("_RipollSanchez2023_midRange", "")
                    t = ax.annotate(
                        label_text,
                        xy=(x_pos + jitter[idx], row[col]),
                        fontsize=7,
                        color=type_color[ctype],
                        va="center",
                        zorder=5,
                    )
                    texts.append(t)

    if annotate and texts:
        adjust_text(
            texts, ax=ax,
            force_points=(1.5, 1.5),
            force_text=(0.8, 0.8),
            expand_points=(2.0, 2.0),
            arrowprops=dict(arrowstyle="-", color="0.4", lw=0.6),
        )

    ax.set_xticks(range(len(TYPE_ORDER)))
    ax.set_xticklabels(TYPE_LABELS, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=12)
    if log_scale:
        ax.set_yscale("log")
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(ylim_bottom if ylim_bottom is not None else ymin * 0.7, ymax * 1.5)
        ymin, ymax = ax.get_ylim()
        base_ticks = [1e-6, 1e-5, 5e-5, 0.0001, 0.0002, 0.0005, 0.001, 0.002,
                      0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10]
        major = [v for v in base_ticks if ymin * 0.8 <= v <= ymax * 1.2]
        ax.yaxis.set_major_locator(ticker.FixedLocator(major))
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(
            lambda v, _: f"{v:g}"
        ))
        ax.yaxis.set_minor_locator(ticker.NullLocator())
    else:
        ymin, ymax = ax.get_ylim()
        margin = (ymax - ymin) * 0.08
        ax.set_ylim(ymin - margin, ymax + margin)
    ax.grid(axis="y", alpha=0.3, linewidth=0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    out_path = f"{OUT_DIR}/{filename}.png"
    fig.savefig(out_path, dpi=100, 
                bbox_inches="tight"
                )
    print(f"Saved {out_path}")
    plt.close(fig)


ANNOT_KWARGS = {
    "density": {"top_n": 3, "bot_n": 0},
}

for col, ylabel, log_scale, filename in METRICS:
    kw = ANNOT_KWARGS.get(filename, {})
    make_strip_plot(col, ylabel, log_scale, filename, annotate=True, **kw)
    make_strip_plot(col, ylabel, log_scale, f"{filename}_clean", annotate=False)


# ── Figure 5: per-connectome heatmap of z-scored metrics ─────────────────────
HEAT_COLS = [
    "density", "mean_weight", "Q", "n_communities",
    "in_degree_cv", "out_degree_cv", "total_degree_cv",
    "reciprocity", "n_reciprocal_pairs", "fraction_in_nontrivial_scc",
    "clustering_coefficient",
]

heat = df.set_index("label")[HEAT_COLS].copy()

# Z-score each column for comparable color scale
heat_z = heat.apply(lambda c: (c - c.mean()) / c.std() if c.std() > 0 else 0, axis=0)

fig2, ax2 = plt.subplots(figsize=(10, 14))
im = ax2.imshow(heat_z.values, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)

ax2.set_yticks(range(len(heat_z)))
ax2.set_yticklabels(heat_z.index, fontsize=6)
ax2.set_xticks(range(len(HEAT_COLS)))
ax2.set_xticklabels(HEAT_COLS, rotation=45, ha="right", fontsize=8)

# Color-coded type bars on the left
for i, label in enumerate(heat_z.index):
    ctype = df.loc[df["label"] == label, "connectome_type"].values[0]
    ax2.plot(-0.8, i, "s", color=type_color[ctype], markersize=5, clip_on=False)

cb = fig2.colorbar(im, ax=ax2, shrink=0.5, label="Z-score")
ax2.set_title("Z-scored Network Metrics per Connectome", fontsize=12, pad=10)
fig2.tight_layout()
fig2.savefig(f"{OUT_DIR}/preassembled_metrics_heatmap.png", dpi=200, bbox_inches="tight")
print(f"Saved {OUT_DIR}/preassembled_metrics_heatmap.png")

plt.close("all")
