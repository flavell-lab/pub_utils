# Plan: Structural Constraints on Molecular Connectomes

## Motivation

Molecular signaling (neurotransmitters, neuropeptides) can be constrained by physical proximity. A neuron can only receive a signal if it's structurally close to the releasing neuron. This plan outlines how to apply structural connectivity as a constraint on molecular connectomes.

## Core Logic

```
constrained[source, target] = molecular[source, target] AND structural[source, target]
```

Where:
- `molecular[source, target]` = 1 if source releases molecule AND target has receptor
- `structural[source, target]` = 1 if source is structurally connected to target

## Design Considerations

### 1. Which structural connectomes?

Available in `connectomes/preassembled/`:

| Type | Datasets | Notes |
|------|----------|-------|
| Chemical synapses | Cook2019, Varshney2011, White*, Witvliet*, Yim2024 | Directional |
| Electrical (gap junctions) | Cook2019, Varshney2011, White*, Witvliet* | Bidirectional |

**Decision needed:** Should the constraint use:
- Chemical only (directional signaling)
- Electrical only (local diffusion)
- Union of both (any structural contact)
- User's choice (parameter)

### 2. Constraint modes

**Mode A: Direct connection (1-hop)**
```
constrained = molecular * (structural > 0)
```
Only allow molecular connections where there's a direct structural synapse.

**Mode B: Multi-hop (n-hop)**
```
reachable = compute_reachability(structural, max_hops=n)
constrained = molecular * reachable
```
Allow molecular connections if target is within n structural hops of source.

**Mode C: Weighted by structural strength**
```
constrained = molecular * structural  # preserves synapse counts
```
Weight molecular connections by structural synapse count.

### 3. API Options

**Option A: Post-processing function**
```python
def constrain_by_structure(
    molecular: pd.DataFrame,
    structural: pd.DataFrame,
    mode: str = 'binary'  # 'binary', 'weighted', 'n-hop'
) -> pd.DataFrame:
    ...
```

**Option B: Integrated into assembly**
```python
assemble_nt_connectome(
    ...,
    structural_constraint='chemical:Cook2019',  # or None for no constraint
    constraint_mode='binary'
)
```

**Option C: Structural connectome loader + constraint function**
```python
structural = load_structural_connectome('chemical', 'Cook2019')
molecular = assemble_nt_connectome(...)
constrained = apply_structural_constraint(molecular, structural)
```

### 4. Handling neuron mismatches

Structural connectomes may have different neuron sets (e.g., male-specific neurons, missing neurons). Need to:
- Align indices before applying constraint
- Warn about neurons present in one but not the other
- Decision: treat missing structural data as "no connection" or "unknown"?

## Proposed Implementation

### Phase 1: Load structural connectomes

```python
def load_structural_connectome(
    synapse_type: str,  # 'chemical', 'electrical', 'both'
    dataset: str,       # 'Cook2019', 'Varshney2011', etc.
    neuron_order: list = AllHermNeurons
) -> pd.DataFrame:
    """Load and standardize a structural connectome."""
    ...
```

### Phase 2: Apply constraint

```python
def apply_structural_constraint(
    molecular: pd.DataFrame,
    structural: pd.DataFrame,
    mode: str = 'binary',
    max_hops: int = 1
) -> pd.DataFrame:
    """
    Constrain molecular connectome by structural connectivity.

    Args:
        molecular: Molecular connectome (source × target)
        structural: Structural connectome (source × target)
        mode: 'binary' (mask), 'weighted' (multiply), 'reachability' (n-hop)
        max_hops: For 'reachability' mode, maximum path length

    Returns:
        Constrained molecular connectome
    """
    ...
```

### Phase 3: Reachability computation (for multi-hop)

```python
def compute_reachability(
    structural: pd.DataFrame,
    max_hops: int
) -> pd.DataFrame:
    """
    Compute which neurons are reachable within n hops.

    Uses matrix exponentiation: reachable = sum(A^k for k in 1..n) > 0
    """
    ...
```

## Metadata extension

When saving constrained connectomes, add to metadata:
```json
{
    "structural_constraint": {
        "dataset": "Cook2019",
        "synapse_type": "chemical",
        "mode": "binary",
        "max_hops": 1
    }
}
```

## Open Questions

1. **Default structural dataset?** Cook2019 is most complete, but should we require explicit choice?

2. **How to handle electrical synapses?** They're bidirectional - should we symmetrize them before applying constraint?

3. **Multi-hop performance?** Matrix exponentiation can be slow for large matrices. Cache reachability matrices?

4. **Biological interpretation?** What hop count corresponds to short/mid/long range in RipollSanchez2023?

## Validation

Compare constrained assembly to preassembled RipollSanchez2023 short-range connectome, which uses structural proximity constraints internally.
