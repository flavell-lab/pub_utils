## Assets Tree

```
assets
├── neuron_features
│   └── neuroanatomy
│
├── structural_connectomes
│   └── preassembled
│       ├── electrical_synapse    [14 sources]
│       └── chemical_synapse      [16 sources]
│
├── molecular_connectomes
│   ├── preassembled
│   │   ├── neuropeptide          [3 sources]
│   │   └── monoamine             [5 sources]
│   │
│   └── candy_assembly (customized logic)
│       ├── neuropeptide
│       ├── classical neurotransmitters
│       └── monoamine
│
├── pairing_info
│   ├── neurotransmitter
│   └── neuropeptide              [3 sources]
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
    │   │   └── reporter
    │   ├── gaba
    │   │   ├── sequencing
    │   │   └── reporter
    │   ├── glutamate
    │   │   └── sequencing
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

## Data Directory Tree

```
data/
├── assets.json
│
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


## External links 

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


## Rationales behind Extract, Transform, Load procedures for molecular connectomes
1) Built monoamine (MA) connectomes and neurotransmitter (NT) connectomes based on fluorescent reporter data for both ligand and receptors
2) Used neuropeptide (NPP) connectomes built by LipollSanchez2023, which was based on in vitro validation and scRNAseq
3) For all connectomes, `1` represents known connection, whereas `0` is the absence of evidence of connection -- not evidence of absence!

### Overlapping information

Data from the same lab - only the latest version was used
e.g. Neuropeptide connectome from LipollSanchez2023 is used instead of Bentley2016_PEP since both came out of the Schaefer lab
e.g. Neurotransmitter atlas from Wang2024 is used instead of Serrano-Saiz2013, Pereira2015, Gendrel2016 since they all came out of the Hobert lab

Data from different labs - all were included as independent observations
exception#1: White1986 (N2U) was dropped and Varshney2011 was used instead because Varshney went back to the lab notebook of White et al. and fixed a few misannotations/ missing links. Varshley also merged N2U with the midbody and tail datasets from 2 other labs.

exception#2: White1986 (A, which was a merge of N2U, N2T, N2W, JSA, JSE)/ Varshley2011 were considered independent observations from Cook2019Herm even though it was the same 5 specimen. This is because Cook2019Herm was a re-annotation of the same images on a new graphical user interface that boosted synaptic counts.