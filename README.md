## Repository organization
```
pub_utils/
├── data/
│   ├── Altun2013/
│   │   └── NT_receptor_info.csv
│   └── RipollSanchez2023/
│       ├── GPCR_per_neuron.csv
│       ├── monoamine_connectome_08062023.csv
│       ├── neuroanatomy.csv
│       ├── NPP_connectome_long_range_01022024.csv
│       ├── NPP_connectome_mid_range_01022024.csv
│       ├── NPP_connectome_short_range_01022024.csv
│       ├── NPP_GPCR_info.csv
│       ├── NPPpairsbyneuron_longrange.csv
│       ├── NPPpairsbyneuron_midrange.csv
│       ├── NPPpairsbyneuron_shortrange.csv
│       └── NPP_per_neuron.csv
├── notebook/
│   ├── analyze_molecular_connectomes.ipynb
│   ├── extract_transform_load_data.ipynb
│   ├── plot_molecular_connectomes.ipynb
│   ├── plot_neuron_features.ipynb
│   └── plot_NPP_GPCR_pairing.ipynb
├── plots/
│   ├── monoamine_connectome.png
│   ├── neuron_features.png
│   ├── NPP_connectome_long_range.png
│   ├── NPP_connectome_mid_range.png
│   ├── NPP_connectome_short_range.png
│   ├── NPP_GPCR_edge_counts.png
│   └── NPP_GPCR_pairing.png
├── processed/
│   ├── monoamine_connectome.pkl
│   ├── neuron_features.pkl
│   ├── NPP_connectome_long_range.pkl
│   ├── NPP_connectome_mid_range.pkl
│   └── NPP_connectome_short_range.pkl
└── src/
    └── pub_utils/
        ├── __init__.py
        ├── core.py
        ├── io.py
        └── plot.py
```

## Data Sources

#### [Cell type classification]
Worm Wiring project S.Cook, Nature 2019 paper
Worm Atlas

#### [Neuroanatomy]
D.Witvliet 2020 bioarxiv paper

#### [Neuropeptide receptor expression]
Worm Atlas - ZF Altun annotation 2013: https://www.wormatlas.org/NTRmainframe.htm
Ripoll-Sanchez...Schaeffer, 2023: https://github.com/LidiaRipollSanchez/NemaMod/tree/main, https://github.com/LidiaRipollSanchez/Neuropeptide-Connectome

#### [Neurotransmitter ligand expression]
Wang...Hobert, 2025: https://pmc.ncbi.nlm.nih.gov/articles/PMC11488851/#s6

#### [Neurotransmitter receptor expression]
Worm Atlas - ZF Altun annotation 2013: https://www.wormatlas.org/NTRmainframe.htm
(pending scRNA data integration)

### [Structural connectome]
OpenWorm C. elegans Connectome Toolbox: https://openworm.org/ConnectomeToolbox/
