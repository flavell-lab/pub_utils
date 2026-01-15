# Molecular Connectome Assembly Plan

## Overview

This document details the implementation plan for assembling molecular connectomes from release and receptor expression data. The assembly logic combines neurotransmitter/neuropeptide release markers with receptor expression to infer functional connectivity.

## Core Assembly Logic

```
Connectome[source, target] = Release[source, NT] AND Receptor[target, receptor]
where (NT, receptor) is a valid pair from pairing_info with confidence >= 1
```

## Data Flow

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                           INPUT DATA LAYER                                   │
├─────────────────────────────────────────────────────────────────────────────┤
│  gene_info.csv          pairing_info.csv       release/*.csv   receptor/*.csv│
│  (NT → genes)           (NT ↔ receptors)       (neuron × gene) (neuron × rec)│
└────────┬────────────────────────┬──────────────────┬───────────────┬────────┘
         │                        │                  │               │
         ▼                        ▼                  ▼               ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                         HELPER FUNCTIONS                                     │
├─────────────────────────────────────────────────────────────────────────────┤
│  _load_assets()         _get_valid_pairs()    _load_release()  _load_receptor()
│  _resolve_markers()     _apply_confidence()   _apply_and_gate() _apply_gate()
└────────┬────────────────────────┬──────────────────┬───────────────┬────────┘
         │                        │                  │               │
         ▼                        ▼                  ▼               ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                         PUBLIC API FUNCTIONS                                 │
├─────────────────────────────────────────────────────────────────────────────┤
│  get_release_vector()                  get_receptor_matrix()                 │
│  → pd.Series (neuron → 0/1/NaN)        → pd.DataFrame (neuron × receptor)    │
└────────────────────────┬───────────────────────────────┬────────────────────┘
                         │                               │
                         ▼                               ▼
┌─────────────────────────────────────────────────────────────────────────────┐
│                       ASSEMBLY FUNCTION                                      │
├─────────────────────────────────────────────────────────────────────────────┤
│  assemble_nt_connectome() / assemble_npp_connectome()                       │
│  → 'per_pair': dict[receptor, DataFrame]                                    │
│  → 'count': DataFrame with receptor pair counts                             │
│  → 'binary': DataFrame with 0/1 values                                      │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## Implementation Steps

### Step 1: Data Loading Infrastructure

#### 1.1 Asset Path Resolution

```python
# File: src/pub_utils/assemble.py

import json
from pathlib import Path

_ASSETS_CACHE = None
_BASE_PATH = Path(__file__).parent.parent.parent.parent  # points to repo root

def _load_assets() -> dict:
    """Load and cache assets.json."""
    global _ASSETS_CACHE
    if _ASSETS_CACHE is None:
        with open(_BASE_PATH / "data" / "assets.json") as f:
            _ASSETS_CACHE = json.load(f)
    return _ASSETS_CACHE

def _get_path(asset_path: str) -> Path:
    """Convert relative asset path to absolute path."""
    return _BASE_PATH / asset_path
```

#### 1.2 Gene Info Loader

```python
_GENE_INFO_CACHE = None

def _load_gene_info() -> pd.DataFrame:
    """Load NT_uptake_synthesis_release_gene_info.csv with lowercase NT names."""
    global _GENE_INFO_CACHE
    if _GENE_INFO_CACHE is None:
        assets = _load_assets()
        path = _get_path(assets["gene_info"]["neurotransmitter"])
        df = pd.read_csv(path)
        df["Neurotransmitter"] = df["Neurotransmitter"].str.lower()
        _GENE_INFO_CACHE = df
    return _GENE_INFO_CACHE
```

#### 1.3 Pairing Info Loader

```python
_PAIRING_INFO_CACHE = {}

def _load_pairing_info(molecule_type: str = "neurotransmitter", source: str = None) -> pd.DataFrame:
    """
    Load receptor-ligand pairing info.

    Args:
        molecule_type: 'neurotransmitter' or 'neuropeptide'
        source: For NPP, one of 'Altun2013', 'Bentley2016', 'RipollSanchez2023'
    """
    cache_key = (molecule_type, source)
    if cache_key not in _PAIRING_INFO_CACHE:
        assets = _load_assets()
        if molecule_type == "neurotransmitter":
            path = _get_path(assets["pairing_info"]["neurotransmitter"])
        else:
            path = _get_path(assets["pairing_info"]["neuropeptide"][source])
        df = pd.read_csv(path)
        df["ligand"] = df["ligand"].str.lower()
        _PAIRING_INFO_CACHE[cache_key] = df
    return _PAIRING_INFO_CACHE[cache_key]
```

---

### Step 2: Marker Resolution

```python
def _resolve_markers(
    neurotransmitter: str,
    markers: list[str]
) -> list[str]:
    """
    Resolve functional category names to gene names.

    Examples:
        _resolve_markers('dopamine', ['release']) -> ['cat-1']
        _resolve_markers('dopamine', ['synthesis', 'release']) -> ['cat-2', 'cat-1']
        _resolve_markers('dopamine', ['cat-2']) -> ['cat-2']  # pass-through

    Raises:
        ValueError: If NT not found or functional category has no gene defined.
    """
    gene_info = _load_gene_info()
    nt_row = gene_info[gene_info["Neurotransmitter"] == neurotransmitter.lower()]

    if nt_row.empty:
        raise ValueError(f"Unknown neurotransmitter: {neurotransmitter}")

    nt_row = nt_row.iloc[0]
    resolved = []

    category_map = {
        "uptake": "UptakeGene",
        "synthesis": "SynthesisGene",
        "release": "ReleaseGene"
    }

    for marker in markers:
        marker_lower = marker.lower()
        if marker_lower in category_map:
            col = category_map[marker_lower]
            gene = nt_row[col]
            if pd.isna(gene) or gene == "":
                raise ValueError(
                    f"No {marker} gene defined for {neurotransmitter}"
                )
            resolved.append(gene)
        else:
            # Assume it's a direct gene name
            resolved.append(marker)

    return resolved
```

---

### Step 3: Release Data Loading

```python
def _load_release_data(
    molecule_type: str,  # 'neurotransmitter' or 'neuropeptide'
    source_method: str,  # 'literature', 'reporter', 'staining', 'sequencing'
    source_dataset: str  # 'Bentley2016', 'Wang2024', 'RipollSanchez2023'
) -> pd.DataFrame:
    """
    Load release expression matrix (neuronID × gene/neuropeptide).

    Returns DataFrame with neuronID as index, genes as columns, values 0/1.
    """
    assets = _load_assets()

    try:
        path_str = assets["release"][molecule_type][source_method][source_dataset]
    except KeyError:
        raise ValueError(
            f"No release data for {molecule_type}/{source_method}/{source_dataset}"
        )

    df = pd.read_csv(_get_path(path_str))
    df = df.set_index("neuronID")
    return df
```

---

### Step 4: Implement `get_release_vector()`

```python
def get_release_vector(
    neurotransmitter: str,
    markers: list[str],
    sources: list[str] = None,
    neuron_order: list = AllHermNeurons
) -> pd.Series:
    """
    Returns neuron vector with release status for a neurotransmitter.

    Args:
        neurotransmitter: 'acetylcholine', 'dopamine', etc.
        markers: ['release'], ['synthesis', 'release'], or specific genes ['cat-2']
        sources: List of source keys like ['literature:Bentley2016', 'reporter:Wang2024']
                 If None, uses all available sources for the NT.
        neuron_order: Canonical neuron ordering for output.

    Returns:
        pd.Series indexed by neuron_order.
        Values: 1 if ALL markers positive (AND gate), 0 if any marker 0, NaN if data missing.

    Logic:
        1. Resolve functional marker names to gene names
        2. Load release matrices from specified sources
        3. For each source, extract columns matching resolved genes
        4. Apply AND gate across genes within each source
        5. Apply OR gate across sources (any positive source = positive)
        6. Align to neuron_order, fill missing with NaN
    """
    # Step 1: Resolve markers
    genes = _resolve_markers(neurotransmitter, markers)

    # Step 2: Determine sources
    if sources is None:
        sources = _get_available_release_sources("neurotransmitter")

    # Step 3: Load and combine data
    combined = pd.DataFrame(index=neuron_order)

    for source_key in sources:
        method, dataset = source_key.split(":")
        try:
            df = _load_release_data("neurotransmitter", method, dataset)
        except ValueError:
            continue

        # Check which genes are in this dataset
        available_genes = [g for g in genes if g in df.columns]
        if not available_genes:
            continue

        # AND gate across genes for this source
        source_result = df[available_genes].min(axis=1)  # min = AND for 0/1
        combined[source_key] = source_result

    # Step 4: OR gate across sources
    if combined.empty:
        return pd.Series(index=neuron_order, dtype=float)

    # max = OR for 0/1, but we need to handle NaN carefully
    # Result is 1 if any source is 1, 0 if all sources are 0, NaN if all NaN
    result = combined.max(axis=1)

    # Align to neuron_order
    result = result.reindex(neuron_order)

    return result
```

---

### Step 5: Receptor Data Loading

```python
def _load_receptor_data(
    molecule_type: str,
    neurotransmitter: str,  # For NT-specific receptor files
    receptor_type: str,     # 'ionotropic', 'metabotropic', 'all'
    source_method: str,     # 'sequencing', 'reporter', 'literature'
    source_dataset: str     # 'Fenyves2020', 'HobertLab', 'Bentley2016', etc.
) -> pd.DataFrame:
    """
    Load receptor expression matrix (neuronID × receptor).

    Navigation through assets.json receptor structure:
        receptor[molecule_type][NT][receptor_type][source_method][source_dataset]
    """
    assets = _load_assets()

    # Navigate the nested structure
    receptor_assets = assets["receptor"][molecule_type]

    # Handle "all" NT case (e.g., Bentley2016 has all NTs in one file)
    nt_key = neurotransmitter.lower()
    if nt_key not in receptor_assets and "all" in receptor_assets:
        nt_key = "all"

    try:
        nt_assets = receptor_assets[nt_key]
        type_assets = nt_assets[receptor_type]
        method_assets = type_assets[source_method]
        path_str = method_assets[source_dataset]
    except KeyError:
        raise ValueError(
            f"No receptor data for {molecule_type}/{nt_key}/{receptor_type}/{source_method}/{source_dataset}"
        )

    df = pd.read_csv(_get_path(path_str))
    df = df.set_index("neuronID")
    return df
```

---

### Step 6: Implement `get_receptor_matrix()`

```python
def get_receptor_matrix(
    neurotransmitter: str,
    sources: list[str],
    gate: str = 'or',
    receptor_type: str = 'all',
    neuron_order: list = AllHermNeurons
) -> pd.DataFrame:
    """
    Returns neuron × receptor matrix for receptors of specified neurotransmitter.

    Args:
        neurotransmitter: Filter receptors by their ligand (e.g., 'dopamine')
        sources: List of source keys ['sequencing:Fenyves2020', 'reporter:HobertLab']
        gate: 'and' (require all sources) or 'or' (any source sufficient)
        receptor_type: 'all', 'ionotropic', or 'metabotropic'
        neuron_order: Canonical neuron ordering for rows.

    Returns:
        pd.DataFrame with neuron_order as index, receptor names as columns.
        Values: 1/0/NaN

    Logic:
        1. Get valid receptors from pairing_info (confidence >= 1)
        2. Filter by receptor_type if specified
        3. Load receptor matrices from each source
        4. For each receptor, apply gate across sources
        5. Return combined matrix aligned to neuron_order
    """
    # Step 1: Get valid receptors for this NT
    pairing = _load_pairing_info("neurotransmitter")
    nt_pairs = pairing[pairing["ligand"] == neurotransmitter.lower()]

    # Filter by confidence
    nt_pairs = nt_pairs[nt_pairs["confidence"] >= 1]

    # Filter by receptor type
    if receptor_type == "ionotropic":
        nt_pairs = nt_pairs[nt_pairs["channel"] == 1]
    elif receptor_type == "metabotropic":
        nt_pairs = nt_pairs[nt_pairs["GPCR"] == 1]

    valid_receptors = set(nt_pairs["receptor"].unique())

    if not valid_receptors:
        return pd.DataFrame(index=neuron_order)

    # Step 2: Load receptor data from each source
    receptor_data_by_source = {}

    for source_key in sources:
        method, dataset = source_key.split(":")

        # Try loading with specific receptor_type, fall back to 'all'
        for rtype in [receptor_type, 'all']:
            try:
                df = _load_receptor_data(
                    "neurotransmitter", neurotransmitter, rtype, method, dataset
                )
                # Filter to valid receptors
                available = [r for r in valid_receptors if r in df.columns]
                if available:
                    receptor_data_by_source[source_key] = df[available]
                    break
            except ValueError:
                continue

    if not receptor_data_by_source:
        return pd.DataFrame(index=neuron_order)

    # Step 3: Combine across sources per receptor
    all_receptors = set()
    for df in receptor_data_by_source.values():
        all_receptors.update(df.columns)

    result = pd.DataFrame(index=neuron_order, columns=sorted(all_receptors))

    for receptor in all_receptors:
        receptor_values = []
        for source_key, df in receptor_data_by_source.items():
            if receptor in df.columns:
                receptor_values.append(df[receptor])

        if not receptor_values:
            continue

        combined = pd.concat(receptor_values, axis=1)

        if gate == 'or':
            # max = OR: 1 if any source is 1
            result[receptor] = combined.max(axis=1).reindex(neuron_order)
        else:  # 'and'
            # min = AND: 1 only if all sources are 1
            result[receptor] = combined.min(axis=1).reindex(neuron_order)

    return result
```

---

### Step 7: Implement `assemble_nt_connectome()`

```python
def assemble_nt_connectome(
    neurotransmitter: str,
    release_markers: list[str],
    release_sources: list[str] = None,
    receptor_sources: list[str] = ['sequencing:Fenyves2020'],
    receptor_gate: str = 'or',
    receptor_type: str = 'all',
    output_format: str = 'binary',
    neuron_order: list = AllHermNeurons
) -> pd.DataFrame | dict[str, pd.DataFrame]:
    """
    Assemble a molecular connectome for a neurotransmitter.

    Args:
        neurotransmitter: e.g., 'dopamine', 'acetylcholine'
        release_markers: ['release'], ['synthesis', 'release'], etc.
        release_sources: Sources for release data, None = all available
        receptor_sources: Sources for receptor data
        receptor_gate: 'or' or 'and' across receptor sources
        receptor_type: 'all', 'ionotropic', 'metabotropic'
        output_format:
            'per_pair' - dict of {receptor: source×target DataFrame}
            'count' - source×target DataFrame with receptor counts
            'binary' - source×target DataFrame with 1 if any connection
        neuron_order: Canonical neuron ordering

    Returns:
        DataFrame or dict depending on output_format.
        Matrix semantics: rows = source neurons, columns = target neurons.
    """
    # Step 1: Get release vector (source capability)
    release = get_release_vector(
        neurotransmitter, release_markers, release_sources, neuron_order
    )

    # Step 2: Get receptor matrix (target capability)
    receptor = get_receptor_matrix(
        neurotransmitter, receptor_sources, receptor_gate, receptor_type, neuron_order
    )

    if receptor.empty:
        if output_format == 'per_pair':
            return {}
        else:
            return pd.DataFrame(
                0, index=neuron_order, columns=neuron_order, dtype=float
            )

    # Step 3: Compute per-receptor connectomes
    per_pair = {}

    for receptor_name in receptor.columns:
        receptor_vec = receptor[receptor_name]

        # Outer product: release[source] × receptor[target]
        # Connection exists if source releases AND target has receptor
        # Using numpy broadcasting: (n,1) * (1,n) = (n,n)
        release_arr = release.values.reshape(-1, 1)
        receptor_arr = receptor_vec.values.reshape(1, -1)

        # Multiply handles NaN propagation automatically
        conn_matrix = release_arr * receptor_arr

        per_pair[receptor_name] = pd.DataFrame(
            conn_matrix,
            index=neuron_order,
            columns=neuron_order
        )

    # Step 4: Format output
    if output_format == 'per_pair':
        return per_pair

    elif output_format == 'count':
        # Sum across receptors, treating NaN as 0 for summation
        count_matrix = pd.DataFrame(
            0.0, index=neuron_order, columns=neuron_order
        )
        for receptor_name, matrix in per_pair.items():
            count_matrix = count_matrix.add(matrix.fillna(0))

        # Restore NaN where ALL receptor matrices were NaN
        all_nan_mask = pd.DataFrame(True, index=neuron_order, columns=neuron_order)
        for matrix in per_pair.values():
            all_nan_mask = all_nan_mask & matrix.isna()
        count_matrix[all_nan_mask] = np.nan

        return count_matrix

    else:  # 'binary'
        count_matrix = assemble_nt_connectome(
            neurotransmitter, release_markers, release_sources,
            receptor_sources, receptor_gate, receptor_type,
            'count', neuron_order
        )
        binary_matrix = (count_matrix >= 1).astype(float)
        binary_matrix[count_matrix.isna()] = np.nan
        return binary_matrix
```

---

### Step 8: Neuropeptide Assembly Function

```python
def assemble_npp_connectome(
    neuropeptide: str = None,  # Specific NPP or None for all
    release_sources: list[str] = ['sequencing:RipollSanchez2023'],
    receptor_sources: list[str] = ['sequencing:RipollSanchez2023'],
    receptor_gate: str = 'or',
    pairing_source: str = 'RipollSanchez2023',
    output_format: str = 'binary',
    neuron_order: list = AllHermNeurons
) -> pd.DataFrame | dict[str, pd.DataFrame]:
    """
    Assemble a neuropeptide connectome.

    Similar to assemble_nt_connectome but:
    - Uses NPP-specific pairing info
    - Neuropeptides are listed directly as column names in release data
    - Can assemble for a single NPP or aggregate all
    """
    # Implementation follows same pattern as NT assembly
    # Key difference: no marker resolution needed, NPP names are direct columns
    pass  # Similar implementation
```

---

## Helper Function Summary

| Function | Purpose | Input | Output |
|----------|---------|-------|--------|
| `_load_assets()` | Load assets.json | - | dict |
| `_get_path(path)` | Resolve relative to absolute path | str | Path |
| `_load_gene_info()` | Load NT gene mapping | - | DataFrame |
| `_load_pairing_info(type, source)` | Load receptor-ligand pairs | str, str | DataFrame |
| `_resolve_markers(nt, markers)` | Convert categories to genes | str, list | list[str] |
| `_load_release_data(type, method, dataset)` | Load release matrix | str×3 | DataFrame |
| `_load_receptor_data(type, nt, rtype, method, dataset)` | Load receptor matrix | str×5 | DataFrame |
| `_get_available_release_sources(type)` | List available sources | str | list[str] |

---

## Testing Plan

### Unit Tests

1. **Marker Resolution**
   - Test all 8 NTs with 'release', 'synthesis', 'uptake'
   - Test pass-through of direct gene names
   - Test error on missing category (e.g., Glutamate uptake)

2. **Release Vector**
   - Verify AND gate: dopamine needs cat-2 AND cat-1
   - Verify OR across sources
   - Check NaN handling for unstudied neurons

3. **Receptor Matrix**
   - Verify receptor filtering by ligand
   - Verify receptor_type filtering (ionotropic vs metabotropic)
   - Test AND vs OR gate across sources

4. **Assembly**
   - Compare output with preassembled Bentley2016 connectomes
   - Verify matrix dimensions (302 × 302)
   - Check that known dopaminergic neurons are sources

### Integration Tests

```python
def test_dopamine_assembly():
    """Verify dopamine connectome against known biology."""
    conn = assemble_nt_connectome(
        'dopamine',
        release_markers=['synthesis', 'release'],  # cat-2 AND cat-1
        receptor_sources=['sequencing:Muralidhara2025'],
        output_format='binary'
    )

    # Known dopaminergic neurons should be sources
    dopaminergic = ['ADEL', 'ADER', 'CEPDL', 'CEPDR', 'CEPVL', 'CEPVR', 'PDEL', 'PDER']
    for neuron in dopaminergic:
        assert conn.loc[neuron].sum() > 0, f"{neuron} should have outgoing connections"

    # Non-dopaminergic neurons should not be sources
    assert conn.loc['ASEL'].sum() == 0, "ASEL should not release dopamine"
```

---

## Output Storage

Assembled connectomes saved to `molecular_connectomes.dk_assembly` in `assets.json`:

```json
{
  "molecular_connectomes": {
    "dk_assembly": {
      "dopamine_synthesis_release_Muralidhara2025": "connectomes/dk_assembly/dopamine_syn_rel_Muralidhara2025.csv",
      "acetylcholine_release_Fenyves2020_ionotropic": "connectomes/dk_assembly/ach_rel_Fenyves2020_ion.csv"
    }
  }
}
```

Naming convention: `{NT}_{markers}_{receptor_source}_{receptor_type}.csv`

---

## Implementation Order

1. **Phase 1: Infrastructure** (data loading)
   - `_load_assets()`, `_get_path()`
   - `_load_gene_info()`, `_load_pairing_info()`
   - `_load_release_data()`, `_load_receptor_data()`

2. **Phase 2: Core Functions**
   - `_resolve_markers()`
   - `get_release_vector()`
   - `get_receptor_matrix()`

3. **Phase 3: Assembly**
   - `assemble_nt_connectome()` with all output formats
   - `assemble_npp_connectome()`

4. **Phase 4: Validation**
   - Unit tests for each function
   - Integration test comparing to Bentley2016 preassembled
   - Documentation and examples

---

## Open Questions for User

1. **Source specification format**: Currently using `"method:dataset"` strings (e.g., `"sequencing:Fenyves2020"`). Is this intuitive, or prefer separate arguments?

2. **Default sources**: Should there be sensible defaults for each NT, or require explicit specification?

3. **NaN semantics**: Current plan treats NaN as "unknown". Should we provide option to treat as 0 (conservative) or 1 (liberal)?

4. **Confidence threshold**: Currently hardcoded at >= 1. Should this be a parameter?

5. **Caching**: Should assembled connectomes be automatically cached, or leave caching to user?
