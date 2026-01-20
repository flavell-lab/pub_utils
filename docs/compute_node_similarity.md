# Computing Neuron Similarity in C. elegans

This document describes how to compute pairwise similarity between the 302 hermaphrodite neurons across four dimensions using the `pub_utils.similarity` module.

## Overview

The similarity module computes how alike any two neurons are based on:

1. **Sensory Type Similarity** - Which sensory modalities neurons respond to (13 types)
2. **Anatomical Process Similarity** - Where neurons have physical processes (28 locations)
3. **NT/NPP Release Similarity** - Which neurotransmitters/neuropeptides neurons release
4. **NT/NPP Receptor Similarity** - Which receptors neurons express

All outputs are 302×302 symmetric matrices with values in [0, 1].

## Quick Start

```python
import pub_utils as pu

# Compute sensory type similarity
sensory = pu.compute_sensory_similarity()
print(sensory['similarity'].loc['ASEL', 'ASER'])  # Compare two neurons

# Find neurons similar to AVAL
similar = pu.find_similar_neurons(sensory['similarity'], 'AVAL', top_n=10)
print(similar)

# Compute all similarity dimensions at once
all_sims = pu.compute_all_similarities(save_path='connectomes/similarity/')
```

## Similarity Metrics

### Jaccard Similarity (default for binary features)

```
Jaccard(A, B) = |A ∩ B| / |A ∪ B|
```

- Returns 1.0 when neurons have identical feature profiles
- Returns 0.0 when neurons share no features
- Appropriate for sensory type, process, and release profiles

### Cosine Similarity (for continuous features)

```
Cosine(A, B) = (A · B) / (||A|| × ||B||)
```

- Returns 1.0 when feature vectors point in same direction
- Returns 0.0 when orthogonal
- Appropriate for receptor expression profiles

## API Reference

### Profile Functions

These functions extract feature vectors for each neuron:

```python
# Sensory modality profile (302 × 13)
sensory_profile = pu.get_sensory_profile()

# Anatomical process profile (302 × 28)
process_profile = pu.get_process_profile()

# NT release profile (302 × 8 neurotransmitters)
nt_release = pu.get_nt_release_profile()

# NPP release profile (302 × ~109 neuropeptides)
npp_release = pu.get_npp_release_profile()

# NT receptor profile (302 × ~77 receptors)
nt_receptor = pu.get_nt_receptor_profile()

# NPP receptor profile (302 × ~154 receptors)
npp_receptor = pu.get_npp_receptor_profile()
```

### Similarity Computation

```python
# Sensory similarity
result = pu.compute_sensory_similarity(metric='jaccard')
# Returns: {'similarity': DataFrame, 'profile': DataFrame, 'metadata': dict}

# Anatomical process similarity
result = pu.compute_process_similarity(metric='jaccard')

# Release similarity (NT + NPP combined)
result = pu.compute_release_similarity(molecule_type='all', metric='jaccard')

# Receptor similarity (NT + NPP combined)
result = pu.compute_receptor_similarity(molecule_type='all', metric='jaccard')

# Compute all at once
all_results = pu.compute_all_similarities(save_path='connectomes/similarity/')
```

### Utility Functions

```python
# Find top N similar neurons
similar = pu.find_similar_neurons(
    similarity_matrix=result['similarity'],
    neuron='AVAL',
    top_n=10,
    threshold=0.5  # Optional minimum similarity
)

# Compare similarity across dimensions for a neuron pair
comparison = pu.compare_similarity_dimensions(
    results=all_results,
    neuron_pair=('ASEL', 'ASER')
)
```

## Data Sources

### Sensory Types (13 modalities)

From `connectomes/neuron_features.csv`:
- chemical, odor, pheromone, osmolarity
- oxygen, carbon_dioxide, photon, thermal
- electrical, mechanical, stretch, nociceptive

### Anatomical Processes (28 locations)

From `connectomes/neuron_features.csv`:
- nerve_ring, dorsal_cord, ventral_cord
- amphid_commissure_left/right
- Various sublateral cords and ganglia

### NT Release

From `data/Bentley2016/` and `data/Wang2024/`:
- Acetylcholine, GABA, Glutamate
- Dopamine, Serotonin, Tyramine, Octopamine, Betaine

### NPP Release

From `data/RipollSanchez2023/`:
- ~109 neuropeptide precursors (flp-*, nlp-*, etc.)

### Receptors

From multiple sources:
- Fenyves2020: Ionotropic receptors (ACh, GABA, Glu)
- HobertLab: Metabotropic receptors
- RipollSanchez2023: NPP receptors

## MCP Server

An MCP (Model Context Protocol) server is available for programmatic access:

```bash
# Start the server
python mcp_server/neuron_similarity.py
```

Add to Claude settings:

```json
{
  "mcpServers": {
    "neuron-similarity": {
      "command": "python",
      "args": ["/path/to/pub_utils/mcp_server/neuron_similarity.py"]
    }
  }
}
```

### Available Tools

| Tool | Description |
|------|-------------|
| `get_neuron_features` | Get sensory/process features for a neuron |
| `get_release_profile` | Get NT/NPP release profile |
| `get_receptor_profile` | Get receptor expression profile |
| `compute_similarity` | Compute similarity between two neurons |
| `find_similar` | Find top-N similar neurons |
| `list_available_data` | List available data sources |
| `get_neurons_by_sensory_type` | Get neurons responding to a modality |
| `get_neurons_by_process` | Get neurons with processes at a location |

### Available Resources

| Resource URI | Description |
|--------------|-------------|
| `similarity://sensory` | Full sensory similarity matrix |
| `similarity://process` | Full process similarity matrix |
| `neurons://list` | List of 302 neurons |
| `features://sensory` | Available sensory features |
| `features://process` | Available process features |

## Output Format

Similarity matrices are saved as CSV files:

```
connectomes/similarity/
├── sensory_similarity.csv      # 302×302, Jaccard similarity
├── process_similarity.csv      # 302×302, Jaccard similarity
├── release_similarity.csv      # 302×302, Jaccard similarity
├── receptor_similarity.csv     # 302×302, Jaccard similarity
```

Matrix format:
- Index and columns: AllHermNeurons (302 neuron IDs)
- Values: [0.0, 1.0] where 1.0 = identical profiles
- Diagonal: Always 1.0 (self-similarity)
- Symmetric: `M[i,j] == M[j,i]`

## NaN Handling

- Features with NaN are excluded from similarity computation
- If two neurons have no overlapping valid features, similarity is NaN
- The metadata includes `n_features` for transparency

## Example Workflow

```python
import pub_utils as pu
import pandas as pd

# 1. Compute all similarities
results = pu.compute_all_similarities()

# 2. Find neurons similar to a query
query = 'AVAL'
for dim in ['sensory', 'process', 'release', 'receptor']:
    print(f"\n{dim.upper()} similarity to {query}:")
    similar = pu.find_similar_neurons(results[dim]['similarity'], query, top_n=5)
    print(similar)

# 3. Compare a specific pair across dimensions
pair = ('ASEL', 'ASER')
comparison = pu.compare_similarity_dimensions(results, pair)
print(f"\nSimilarity between {pair[0]} and {pair[1]}:")
print(comparison)

# 4. Save matrices for downstream analysis
for dim, result in results.items():
    result['similarity'].to_csv(f'connectomes/similarity/{dim}.csv')
```

## Reproducibility

The similarity computation is fully deterministic:
- No random components
- Fixed neuron ordering (AllHermNeurons)
- Consistent data loading from assets.json

To regenerate matrices:
```python
pu.compute_all_similarities(save_path='connectomes/similarity/')
```
