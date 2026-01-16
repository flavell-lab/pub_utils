## Project Overview

This repository provides utilities for extracting, transforming, and loading publicly available C. elegans neuroscience datasets.

## Installation (uv)

```bash
git clone https://github.com/flavell-lab/pub_utils.git
cd pub_utils
uv sync
```
The package requires Python >= 3.10. Dependencies are managed in `pyproject.toml`.

## Package Architecture

The `pub_utils` package (`src/pub_utils/`) exports:

- **`NeuronFeatures`** (`core.py`): Stores neuroanatomical features as a matrix (neurons x features). Features are categorized into `cellType`, `sensoryType`, `segment`, and `process`. Provides fast lookup by neuron ID or feature name.

- **`NeuronInteraction`** (`core.py`): Square adjacency matrix wrapper for connectome data. Handles source→recipient relationships, reciprocal pairs detection, BFS shortest path, and degree analysis. Matrix values represent connection strength (e.g., number of unique ligand-receptor pairs).

- **`plot_connectome_matrix`** / **`plot_reciprocal_network`** (`plot.py`): Visualization functions using seaborn heatmaps and networkx graphs with discrete colormaps.

And many other useful functionality -- see `src/pub_utils/__init__.py` for the full list of features.

## Assets Tree

- **Gene Info**: `data/HobertLab/NT_uptake_synthesis_release_gene_info.csv` - maps functional categories (uptake, synthesis, release) to gene names for each NT
- **Pairing Info**: `data/Altun2013/NT_receptor_info.csv` - maps receptors to ligands with confidence scores and receptor type flags (ionotropic/metabotropic)
- **Release Data**: neuron × gene matrices (binary) from literature, reporter, staining methods
- **Receptor Data**: neuron × receptor matrices (binary) from sequencing, reporter, literature methods

File paths in `assets.json`
```
assets
├── neuron_features
│
├── connectomes
│   ├── preassembled/
│   │   ├── structural/           (accessed through OpenWorm)
│   │   │   ├── chemical
│   │   │   └── electrical
│   │   └── molecular/            (from Bentley2016 and RipollSanchez2023)
│   │
│   ├── candy_assembly/           (customized logic, all are molecular)
│   │   ├── dopamine/         
│   │   ├── serotonin/
│   │   ├── tyramine/
│   │   ├── octopamine/
│   │   ├── acetylcholine/
│   │   ├── gaba/
│   │   ├── glutamate/
│   │   ├── individual_neuropeptides/
│   │   ├── aggregated_neuropeptides/
│   │   ├── aggregated_synapticNT/
│   │   └── aggregated_extrasynapticNT/
│   │
│   └── reproducibility_tests/    (validation reports on molecular assembly from claude code)
│             
├── pairing_info
│   ├── neurotransmitter
│   └── neuropeptide             
│
├── release
│   ├── neurotransmitter
│   │   ├── literature
│   │   ├── reporter
│   │   └── staining
│   │
│   └── neuropeptide
│       ├── literature
│       └── sequencing
│
└── receptor
    ├── neurotransmitter
    │   ├── acetylcholine
    │   │   ├── sequencing
    │   │   └── reporter          (metabotropic only)
    │   ├── glutamate
    │   │   ├── sequencing
    │   │   └── reporter          (metabotropic only)
    │   ├── gaba
    │   │   ├── sequencing
    │   │   └── reporter
    │   ├── dopamine
    │   │   ├── reporter
    │   │   └── sequencing
    │   ├── serotonin
    │   │   └── reporter
    │   ├── tyramine
    │   │   └── sequencing
    │   ├── octopamine
    │   │   └── sequencing
    │   └── all
    │       └── literature
    │   
    └── neuropeptide
        ├── literature
        └── sequencing
```
Note: `{username}_assembly/` directories contain custom connectomes assembled via `notebook/assemble_connectomes.ipynb`. Single-molecule neurotransmitter connectomes go into molecule subdirectories; neuropeptide and aggregated connectomes are saved flat in the root.


## Data Directory Tree

```
data/
├── Altun2013/
│   ├── NPP_receptor_info.csv
│   └── NT_receptor_info.csv
│
├── Bentley2016/
│   ├── NPP_receptor_info.csv
│   ├── NPP_receptor_metabotropic_literature.csv
│   ├── NPP_release_literature.csv
│   ├── NT_receptor_all_literature.csv
│   ├── NT_release_literature.csv
│   ├── dopamine_receptor_all_reporter.csv
│   ├── octopamine_receptor_all_literature.csv
│   ├── serotonin_receptor_all_literature.csv
│   ├── tyramine_receptor_all_literature.csv
│   ├── monoamine_expression.csv                   (raw)
│   ├── monoamine_receptor_expression.csv          (raw)
│   ├── neuropeptide_expression.csv                (raw)
│   ├── neuropeptide_receptor_expression.csv       (raw)
│   └── supplementary_references.csv               (raw)
│
├── Dag2023/
│   ├── serotonin_receptor_all_reporter.csv
│   └── 5htr_expression_dv_final.csv               (raw)
│
├── Fenyves2020/
│   ├── acetylcholine_receptor_ionotropic_sequencing.csv
│   ├── gaba_receptor_ionotropic_sequencing.csv
│   ├── glutamate_receptor_ionotropic_sequencing.csv
│   ├── NT_receptor_expression.csv                 (raw)
│   └── NT_receptor_polarity.csv                   (raw)
│
├── HobertLab/
│   ├── NT_uptake_synthesis_release_gene_info.csv
│   ├── acetylcholine_receptor_metabotropic_reporter.csv
│   ├── gaba_receptor_all_reporter.csv
│   ├── MA_gaba_release_expression_sequencing.csv  (raw)
│   └── NT_receptors.R                             (raw)
│
├── Muralidhara2025/
│   ├── dopamine_receptor_all_reporter.csv
│   └── dopamine_receptor_all_sequencing.csv
│
├── RipollSanchez2023/
│   ├── NPP_receptor_info.csv
│   ├── NPP_receptor_all_sequencing.csv
│   ├── NPP_release_sequencing.csv
│   ├── neuroanatomy.csv
│   ├── NPP_connectome_short_range_01022024.csv    (raw)
│   ├── NPP_connectome_mid_range_01022024.csv      (raw)
│   ├── NPP_connectome_long_range_01022024.csv     (raw)
│   ├── monoamine_connectome_08062023.csv          (raw)
│   ├── GPCR_per_neuron.csv                        (raw)
│   ├── NPP_per_neuron.csv                         (raw)
│   ├── NPPpairsbyneuron_*.csv                     (raw)
│   ├── 30072020_CENGEN_*.csv                      (raw)
│   ├── group/                                     (raw)
│   └── individual/                                (raw)
│
└── Wang2024/
    ├── NT_release_reporter.csv
    ├── NT_release_staining.csv
    ├── NT_release_reporter_male.csv
    └── NT_release_staining_male.csv
```
Note: Files marked `(raw)` are original source files kept for reference but not directly used in assets.json.

## Data Conventions

- Connectome matrices are square DataFrames: rows = source neurons, columns = recipient neurons
- Value `1` = known connection; `0.5` = variable connection; `0` = evidence of absence; `NaN` = absence of evidence
- Neuron IDs follow WormAtlas naming (e.g., `ASEL`, `ASER`, `DA01`)

## Asset Construction Workflow

1. **Extract**: Notebooks in `notebook/` convert source formats (RData, XLSX) to CSV
2. **Transform**: One-hot encode categorical features, standardize neuron ordering via `standardize_dataframe(df, neuron_order)`, saved to CSV
3. **Load**: Wrap in `NeuronFeatures` or `NeuronInteraction` classes for analysis

## Molecular Connectome Assembly

### Logic

```
Connectome[source, target] = Release[source, NT] AND Receptor[target, receptor] AND (optional) Constraint
where (NT, receptor) is a valid pair from pairing_info with confidence >= 1
and 
Constraint is set with structural chemical synapses and applicable to classical NT (acetylcholine, glutamate)
```

### Design

1. **Release filtering**: User chooses AND or OR gate across data sources; AND is recommended
2. **Receptor filtering**: User chooses AND or OR gate across data sources; OR is recommended
3. **Confidence threshold**: Applied in code logic (confidence < 1 treated as 0)
4. **Missing values**: Preserved as NaN throughout, all functions robust to NaN
5. **Output location**: Assembled connectomes saved under `connectomes/{username}_assembly` in `assets.json`


## External links 

#### [Neuroanatomy]
Worm Atlas:
https://www.wormatlas.org/neurons/Individual%20Neurons/Neuronframeset.html

Witvliet 2021:
https://github.com/dwitvliet/nature2021/tree/master

RipollSanchez 2023 (subcellular localization):
https://github.com/LidiaRipollSanchez/Neuropeptide-Connectome

#### [Neuropeptide ligand & receptor]
Worm Atlas - Altun 2013: 
https://www.wormatlas.org/NTRmainframe.htm

RipollSanchez...Schaeffer 2023 (fluorescent reporter & scRNAseq): 
https://github.com/LidiaRipollSanchez/NemaMod/tree/main
https://github.com/LidiaRipollSanchez/Neuropeptide-Connectome

#### [Monoamine/Neurotransmitter ligand]
Wang...Hobert, 2025:
https://pmc.ncbi.nlm.nih.gov/articles/PMC11488851/#s6
https://iiif.elifesciences.org/lax:95402%2Felife-95402-fig3-v1.tif/full/,1500/0/default.jpg

WormAtlas:
https://www.wormatlas.org/neurotransmitterstable.htm

#### [Monoamine/Neurotransmitter receptor]
Worm Atlas - Altun 2013: 
https://www.wormatlas.org/NTRmainframe.htm

GABA-A receptors (fluorescent reporter) - Gendrel...Hobert 2016: 
https://elifesciences.org/articles/17686#tbl4

GABA-B receptors (fluorescent reporter) - Yemini...Hobert 2023:
https://pmc.ncbi.nlm.nih.gov/articles/PMC10494711/#SM1

Dopamine receptors (fluorescent reporter & scRNAseq) - Muralidhara & Hardege 2025: 
https://pmc.ncbi.nlm.nih.gov/articles/PMC12539964/table/T4

Serotonin receptors (fluorescent reporter) - Dag...Flavell 2023: 
CSV curated by Ugur Dag for Di Kang to make Figure 7

### [Structural connectome]
White 1986, Varshley 2011, Cook 2019, Cook 2020, Witvliet 2021 - accessed via OpenWorm C. elegans Connectome Toolbox: 
https://openworm.org/ConnectomeToolbox/
