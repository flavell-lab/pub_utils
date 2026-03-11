# Key Insights from Preassembled Connectome Metrics

Computed from `preassembled_metrics.csv` (131 connectomes: 16 chemical synapse, 15 electrical synapse, 94 neuropeptide, 6 monoamine).

---

## 1. Density and Scale

| Type | Density range | Edge count range | Mean weight range |
|------|--------------|-----------------|-------------------|
| Chemical synapse | 0.0016 -- 0.040 | 147 -- 3,805 | 1.2 -- 9.0 |
| Electrical synapse | 0.0006 -- 0.024 | 56 -- 2,185 | 1.1 -- 7.7 |
| Neuropeptide (aggregated) | 0.44 -- 0.64 | 40,254 -- 58,112 | 2.7 -- 4.3 |
| Neuropeptide (individual) | 0.0002 -- 0.12 | 20 -- 11,162 | 1.0 |
| Monoamine | 0.0006 -- 0.032 | 54 -- 2,866 | 1.0 -- 1.4 |

- **Structural connectomes are sparse** (density < 5%), while **aggregated neuropeptide connectomes are dense** (44--64% of all possible edges present).
- Individual neuropeptide ligand-receptor pair connectomes span 3 orders of magnitude in density, from extremely sparse (FLP-17/FRPR-8: 20 edges) to moderately dense (FLP-9/DMSR-7: 11,162 edges).
- All 92 individual neuropeptide connectomes have **uniform weights of 1.0** (binary connectivity), while aggregated neuropeptide and structural connectomes carry variable weights.

---

## 2. Modularity (Community Structure)

| Type | Q range | Typical n_communities |
|------|---------|----------------------|
| Chemical synapse | 0.19 -- 0.59 | 8 -- 283 |
| Electrical synapse | 0.52 -- 0.84 | 14 -- 284 |
| Neuropeptide (aggregated) | 0.02 -- 0.16 | 2 -- 5 |
| Neuropeptide (individual) | 0.00 -- 0.33 | 112 -- 302 |
| Monoamine | 0.00 -- 0.32 | 19 -- 302 |

- **Electrical synapse networks are the most modular** (mean Q = 0.70), consistent with gap junctions forming tight local clusters.
- **Chemical synapse networks show moderate modularity** (mean Q = 0.47), with the best-reconstructed datasets (Cook2019, Cook2019Male) at Q > 0.5.
- **Aggregated neuropeptide networks have near-zero modularity** (Q = 0.02 for Bentley2016, Q = 0.16 for RipollSanchez2023), reflecting their high density -- when most neurons connect to most others, community structure is destroyed.
- **12 individual neuropeptide pairs** are so sparse that Louvain detects 302 communities (all singletons), meaning no neurons cluster together at all. These are ultra-narrow signaling channels (e.g., NLP-23/GNRR-3 with only 22 edges).

---

## 3. Hubness (Degree Distribution Inequality)

| Type | total_degree_cv range | Interpretation |
|------|----------------------|----------------|
| Chemical synapse | 0.93 -- 4.54 | Moderate to high hub structure |
| Electrical synapse | 1.45 -- 4.71 | Moderate to high |
| Neuropeptide (aggregated) | 0.57 -- 0.64 | Very uniform |
| Neuropeptide (individual) | 1.18 -- 6.66 | Extremely variable |
| Monoamine | 2.02 -- 6.50 | Very hub-dominated |

- **Aggregated neuropeptide networks have the most uniform degree distributions** (CV ~ 0.6), meaning all neurons participate roughly equally -- a consequence of high density.
- **Individual neuropeptide and monoamine connectomes are the most hub-dominated**, with some showing total_degree_cv > 6. This reflects the biology: only a few neurons release a given neuropeptide or express its receptor, creating extreme sender/receiver asymmetry.
- **In-degree vs out-degree asymmetry** is striking for molecular connectomes: out_degree_cv reaches 12.2 (FLP-23/DMSR-7, FLP-1/DMSR-7), indicating that a single neuron (or very few) releases the ligand to many targets. In contrast, in_degree_cv is typically 1.5--5, reflecting broader receptor expression.

---

## 4. Reciprocity and Recurrence

| Type | Reciprocity range | Fraction in non-trivial SCC |
|------|------------------|----------------------------|
| Chemical synapse | 0.17 -- 0.54 | 0.06 -- 0.97 |
| Electrical synapse | 0.87 -- 1.00 | 0.07 -- 0.99 |
| Neuropeptide (aggregated) | 0.62 -- 0.65 | 0.95 |
| Neuropeptide (individual) | 0.00 -- 0.34 | 0.00 -- 0.17 |
| Monoamine | 0.00 -- 0.06 | 0.00 -- 0.06 |

- **Electrical synapses are nearly perfectly reciprocal** (14/15 have reciprocity > 0.99), confirming their bidirectional nature. The sole exception is Cook2019Male (r = 0.87), likely due to reconstruction asymmetries in the male-specific neurons.
- **38 out of 92 individual neuropeptide connectomes have zero reciprocity** -- no neuron pair communicates bidirectionally through that ligand-receptor channel. This is expected for ligand-receptor pairs where releasing and receiving neurons don't overlap.
- **Chemical synapse reciprocity varies substantially across datasets** (0.17 for WitvlietL1_01 to 0.54 for Cook2020), reflecting differences in reconstruction completeness and developmental stage.
- **Aggregated neuropeptide networks have high SCC fractions** (~0.95), meaning nearly all 302 neurons participate in directed cycles when all neuropeptide channels are combined. Individual channels are the opposite: most have SCC fractions near zero.
- **Monoamine signaling is almost entirely unidirectional** (reciprocity < 0.07, SCC fraction < 0.06), consistent with a broadcast signaling mode.

---

## 5. Dataset-level Observations

### Structural datasets (chemical + electrical)
- **Cook2019 hermaphrodite** is the most complete wiring diagram: highest chemical synapse density (0.040), most edges (3,671), and largest giant SCC (97% of neurons).
- **Cook2020** appears to be a partial reconstruction: only 147 chemical synapse edges and 56 electrical edges, with extremely high degree CV (> 4.5) and low SCC coverage (6%).
- **Developmental progression** is visible across the Witvliet series: edge count and SCC fraction increase monotonically from L1_01 through A_08, reflecting circuit maturation.

### Molecular datasets
- **FLP-9/DMSR-7** is the densest individual neuropeptide channel (12.3% density, 11,162 edges), suggesting FLP-9 and DMSR-7 are broadly expressed.
- **FLP-9/DMSR-1** is the second-densest (9.5% density, 8,601 edges), with FLP-9 appearing in both -- consistent with FLP-9 being a widely released neuropeptide.
- **Bentley2016 aggregated monoamine** has fewer edges (1,916) than **RipollSanchez2023 aggregated monoamine** (2,866), but Bentley2016 includes separate per-transmitter breakdowns (dopamine, serotonin, tyramine, octopamine). Dopamine alone accounts for 1,152 edges (60% of Bentley2016 total).

---

## 6. Clustering Coefficient (Transitivity)

The global clustering coefficient $C$ measures the fraction of connected triples that close into triangles. High $C$ means that if neuron A connects to both B and C, then B and C are also likely connected.

| Type | CC range | Triangle count range |
|------|---------|---------------------|
| Chemical synapse | 0.14 -- 0.43 | 256 -- 6,853 |
| Electrical synapse | 0.00 -- 0.19 | 0 -- 807 |
| Neuropeptide (aggregated) | 0.70 -- 0.80 | 1.46M -- 3.14M |
| Neuropeptide (individual) | 0.00 -- 0.66 | 0 -- 244,940 |
| Monoamine | 0.00 -- 0.06 | 0 -- 12,723 |

- **Aggregated neuropeptide networks are extremely clustered** (CC = 0.70--0.80). With 44--64% density, most triples naturally close into triangles. Bentley2016 contains over 3.1 million triangles.
- **Chemical synapse networks show moderate clustering** (mean CC = 0.18). This is biologically meaningful: if neuron A synapses onto both B and C, there is roughly an 18% chance B and C also share a synapse, reflecting local circuit motifs.
- **Cook2020 chemical synapse is an outlier** with the highest structural CC (0.43) despite having the fewest edges (147). This is a statistical artefact of its extreme sparsity: the few edges that exist happen to form a higher proportion of triangles, but the absolute count is small (256 triangles).
- **Electrical synapse clustering is lower than chemical** (mean CC = 0.11 vs 0.18), and increases with developmental stage in the Witvliet series (L1_01: 0.024, L1_04: 0.118, A_07: 0.165), suggesting gap junction triangle motifs emerge during maturation.
- **Monoamine networks have near-zero clustering** (mean CC = 0.035). This is consistent with monoamine signaling being broadcast-like: a releasing neuron signals broadly to many targets that do not connect back to each other. Tyramine has zero triangles entirely.
- **32 out of 92 individual neuropeptide connectomes have zero clustering** -- no triangles exist. These correspond to the sparsest ligand-receptor pairs where connectivity is too narrow for three-node loops to form.
- **FLP-9/EGL-6** has the highest individual neuropeptide CC (0.66), despite only 1,697 edges. The FLP-9 ligand reappears here as it did in the density rankings (FLP-9/DMSR-7, FLP-9/DMSR-1), confirming FLP-9 as a broadly connected signaling system whose targets share extensive interconnection.
- **Clustering and density are positively correlated** overall (Pearson r = 0.63), but this relationship reverses within chemical synapses (r = -0.23), where the densest reconstructions (Cook2019) have lower CC than sparse ones (Cook2020). This suggests that denser structural reconstructions reveal more open triples that don't close, while sparse reconstructions preferentially capture the strongest (and most clustered) connections.
