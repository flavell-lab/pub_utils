# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository provides utilities for extracting, transforming, and loading publicly available C. elegans neuroscience datasets. It focuses on building molecular connectomes (neuropeptide, monoamine, neurotransmitter) and neuroanatomical feature matrices for quantitative modeling.

## Setup and Installation

```bash
# Create virtual environment and install
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

The package requires Python >= 3.10. Dependencies are managed in `pyproject.toml`.

## Package Architecture

The `pub_utils` package (`src/pub_utils/`) exports:

- **`NeuronFeatures`** (`core.py`): Stores neuroanatomical features as a matrix (neurons x features). Features are categorized into `cellType`, `sensoryType`, `segment`, and `process`. Provides fast lookup by neuron ID or feature name.

- **`NeuronInteraction`** (`core.py`): Square adjacency matrix wrapper for connectome data. Handles source→recipient relationships, reciprocal pairs detection, BFS shortest path, and degree analysis. Matrix values represent connection strength (e.g., number of unique ligand-receptor pairs).

- **`plot_connectome_matrix`** / **`plot_reciprocal_network`** (`plot.py`): Visualization functions using seaborn heatmaps and networkx graphs with discrete colormaps.

- **`handle_pickle`** / **`get_file_for_pair`** / **`standardize_dataframe`** (`io.py`): I/O utilities for pickle serialization and neuropeptide-receptor pair file lookups.

- **`constants.py`**: Canonical neuron ID lists - `AllHermNeurons` (302), `AllMaleNeurons` (385), `SexSharedNeurons`, `HermSpecificNeurons`, `MaleSpecificNeurons`.

## Data Conventions

- Connectome matrices are square DataFrames: rows = source neurons, columns = recipient neurons
- Value `1` = known connection; `0` = absence of evidence (not evidence of absence)
- Neuron IDs follow WormAtlas naming (e.g., `ASEL`, `ASER`, `DA01`)
- Processed data stored as pickle files in `processed/`

## Workflow

1. **Extract**: Notebooks in `notebook/` convert source formats (RData, XLSX) to CSV
2. **Transform**: One-hot encode categorical features, standardize neuron ordering via `standardize_dataframe(df, neuron_order)`
3. **Load**: Wrap in `NeuronFeatures` or `NeuronInteraction` classes, serialize to pickle

## Data Sources

Primary datasets from RipollSanchez2023 (neuropeptide connectomes), Wang2024/Hobert lab (neurotransmitter atlas), and OpenWorm Connectome Toolbox (structural connectomes). See README.md for full citation details.
