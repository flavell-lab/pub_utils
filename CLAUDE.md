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

## Claude Code Session Guidelines

1. **Daily logging**: Keep track of output files and the logic behind their creation. Save logs under `claudecode/yyyy-mm-dd.py` (e.g., `claudecode/2026-01-14.py`). Prepare the daily log file when I have <3% context left until auto-compact.

2. **Ask, don't guess**: When requirements are ambiguous or in conflict with each other, ask the user for clarification rather than interpolating or guessing.

3. **Token efficiency**: Use simple queries for easy tasks to save tokens. Avoid over-exploring when the task is straightforward.

4. **Clear documentation**: The logic behind the data ETL and connectome assembly needs to be extremely clearly documented for future contributors to reproduce the outputs.

5. **Continuity between sessions**: At the beginning of each claude session, load in the last log file from `claudecode` so that you know where to pick up from.

6. **Error handling**: Minimize try/except patterns in `src/` functions. If a try block is necessary, emit a `warnings.warn()` message on failure rather than silently continuing. This ensures users are aware when something unexpected happens.

## Molecular Connectome Assembly Plan

### Overview

Assemble molecular connectomes by combining neurotransmitter (NT) release data with receptor expression data, using valid NT-receptor pairings from curated pairing info files.

### Data Assets (see `data/assets.json`)

- **Gene Info**: `data/HobertLab/NT_uptake_synthesis_release_gene_info.csv` - maps functional categories (uptake, synthesis, release) to gene names for each NT
- **Pairing Info**: `data/Altun2013/NT_receptor_info.csv` - maps receptors to ligands with confidence scores and receptor type flags (ionotropic/metabotropic)
- **Release Data**: neuron × gene matrices (binary) from literature, reporter, staining methods
- **Receptor Data**: neuron × receptor matrices (binary) from sequencing, reporter, literature methods

### Assembly Logic

```
Connectome[source, target] = Release[source, NT] AND Receptor[target, receptor]
where (NT, receptor) is a valid pair from pairing_info with confidence >= 1
```

### Key Assembly Design Decisions

1. **Release filtering**: AND gate only - all specified markers must be positive
2. **Receptor filtering**: User chooses AND or OR gate across sources
3. **Confidence threshold**: Applied in code logic (confidence < 1 treated as 0), not in data files
4. **Missing values**: Preserved as NaN throughout, all functions robust to NaN
5. **Output location**: Assembled connectomes saved under `molecular_connectomes.dk_assembly` in assets.json
