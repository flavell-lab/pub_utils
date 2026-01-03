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

## Data 

#### [Neuroanatomy]
Worm Atlas:
https://www.wormatlas.org/neurons/Individual%20Neurons/Neuronframeset.html

Witvliet 2021:
https://github.com/dwitvliet/nature2021/tree/master

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
White 1986, Witvliet 2021 - accessed via OpenWorm C. elegans Connectome Toolbox: 
https://openworm.org/ConnectomeToolbox/


# Rationales behind Extract, Transform, Load procedures for molecular connectomes
1) Built monoamine (MA) connectomes and neurotransmitter (NT) connectomes based on fluorescent reporter data for both ligand and receptors
2) Used neuropeptide (NPP) connectomes built by LipollSanchez2023, which was based on in vitro validation and scRNAseq
3) For all connectomes, `1` represents known connection, whereas `0` is the absence of evidence of connection -- not evidence of absence!