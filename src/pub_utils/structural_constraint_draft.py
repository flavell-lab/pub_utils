"""
DRAFT - Structural Constraints on Molecular Connectomes
========================================================

This is a DRAFT module implementing structural constraints on molecular connectomes.
NOT YET ITERATED ON - created 2026-01-15 based on docs/structural_constraint_plan.md

The goal is to constrain molecular signaling by physical proximity: a neuron can only
receive a signal if it's structurally close to the releasing neuron.

Core logic:
    constrained[source, target] = molecular[source, target] AND structural[source, target]

Key design decisions (from plan):
    - Electrical synapses should be symmetrized (they're bidirectional)
    - Default structural constraint = None (no constraint unless specified)
    - Modes: 'binary' (mask), 'weighted' (multiply)
"""

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from .constants import AllHermNeurons

# Base path points to repository root
_BASE_PATH = Path(__file__).parent.parent.parent


def _load_assets() -> dict:
    """Load assets.json (cached in assemble.py, but we load fresh here for draft)."""
    with open(_BASE_PATH / "data" / "assets.json") as f:
        return json.load(f)


def _get_path(asset_path: str) -> Path:
    """Convert relative asset path to absolute path."""
    return _BASE_PATH / asset_path


def _symmetrize_matrix(df: pd.DataFrame) -> pd.DataFrame:
    """
    Symmetrize a matrix by taking element-wise maximum of (A, A^T).

    This is appropriate for electrical synapses (gap junctions) which are
    bidirectional - if A connects to B, B also connects to A.

    Args:
        df: Square DataFrame

    Returns:
        Symmetrized DataFrame where df[i,j] = df[j,i] = max(original[i,j], original[j,i])
    """
    values = df.values
    symmetric = np.maximum(values, values.T)
    return pd.DataFrame(symmetric, index=df.index, columns=df.columns)


def load_structural_connectome(
    synapse_type: str,
    dataset: str,
    neuron_order: list = None,
    symmetrize_electrical: bool = True
) -> pd.DataFrame:
    """
    Load a structural connectome from preassembled files.

    Args:
        synapse_type: 'chemical', 'electrical', or 'both'
        dataset: Dataset name (e.g., 'Cook2019', 'Varshney2011', 'WhiteWhole')
        neuron_order: Canonical neuron ordering for output. Default: AllHermNeurons
        symmetrize_electrical: If True, symmetrize electrical synapse matrices
                               (gap junctions are bidirectional). Default: True

    Returns:
        pd.DataFrame (source x target) with synapse counts.
        For 'both', returns union (sum) of chemical and electrical.

    Example:
        >>> structural = load_structural_connectome('chemical', 'Cook2019')
        >>> structural = load_structural_connectome('both', 'Cook2019')
    """
    if neuron_order is None:
        neuron_order = AllHermNeurons

    assets = _load_assets()
    chem = None
    elec = None

    if synapse_type in ['chemical', 'both']:
        chem_assets = assets['structural_connectomes']['preassembled']['chemical_synapse']
        chem_path = chem_assets.get(dataset)
        if chem_path is None:
            available = list(chem_assets.keys())
            raise ValueError(
                f"No chemical synapse data for dataset '{dataset}'. "
                f"Available: {available}"
            )
        chem = pd.read_csv(_get_path(chem_path), index_col=0)
        chem = chem.reindex(index=neuron_order, columns=neuron_order).fillna(0)

    if synapse_type in ['electrical', 'both']:
        elec_assets = assets['structural_connectomes']['preassembled']['electrical_synapse']
        elec_path = elec_assets.get(dataset)
        if elec_path is None:
            available = list(elec_assets.keys())
            raise ValueError(
                f"No electrical synapse data for dataset '{dataset}'. "
                f"Available: {available}"
            )
        elec = pd.read_csv(_get_path(elec_path), index_col=0)
        elec = elec.reindex(index=neuron_order, columns=neuron_order).fillna(0)

        # Symmetrize electrical synapses (gap junctions are bidirectional)
        if symmetrize_electrical:
            elec = _symmetrize_matrix(elec)

    if synapse_type == 'chemical':
        return chem
    elif synapse_type == 'electrical':
        return elec
    else:  # 'both'
        return chem + elec


def apply_structural_constraint(
    molecular: pd.DataFrame,
    structural: pd.DataFrame,
    mode: str = 'binary'
) -> pd.DataFrame:
    """
    Constrain molecular connectome by structural connectivity.

    Args:
        molecular: Molecular connectome (source x target)
        structural: Structural connectome (source x target),
                    e.g., from load_structural_connectome() or compute_reachability_matrix()
        mode:
            'binary' - Mask: keep molecular connection only if structural connection exists
                       Result = molecular * (structural > 0)
            'weighted' - Multiply: weight molecular by structural synapse count
                         Result = molecular * structural

    Returns:
        Constrained molecular connectome with same shape as input.
        Neurons missing from structural are treated as having no connections (0).

    Example:
        >>> molecular = pu.assemble_nt_connectome('dopamine', ...)['binary']
        >>> structural = load_structural_connectome('chemical', 'Cook2019')
        >>> constrained = apply_structural_constraint(molecular, structural)
    """
    # Align indices
    common_idx = molecular.index.intersection(structural.index)
    common_col = molecular.columns.intersection(structural.columns)

    if len(common_idx) < len(molecular.index):
        missing = set(molecular.index) - set(common_idx)
        warnings.warn(
            f"{len(missing)} neurons in molecular but not structural connectome. "
            f"First 5: {list(missing)[:5]}. These will be treated as having no structural connections."
        )

    mol_aligned = molecular.loc[common_idx, common_col]
    struct_aligned = structural.loc[common_idx, common_col]

    if mode == 'binary':
        # Mask: molecular * (structural > 0)
        mask = (struct_aligned > 0).astype(float)
        constrained = mol_aligned * mask
    elif mode == 'weighted':
        # Weighted: molecular * structural
        constrained = mol_aligned * struct_aligned
    else:
        raise ValueError(f"Unknown mode '{mode}'. Expected 'binary' or 'weighted'.")

    # Restore full shape with zeros for missing neurons
    result = pd.DataFrame(0.0, index=molecular.index, columns=molecular.columns)
    result.loc[common_idx, common_col] = constrained

    return result


def constrain_molecular_connectome(
    molecular: pd.DataFrame,
    structural_dataset: str = 'Cook2019',
    synapse_type: str = 'both',
    mode: str = 'binary',
    symmetrize_electrical: bool = True
) -> dict:
    """
    High-level function to constrain a molecular connectome by structural connectivity.

    This is a convenience wrapper that loads the structural connectome and applies
    the constraint in one step.

    Args:
        molecular: Molecular connectome (source x target), e.g., from assemble_nt_connectome()
        structural_dataset: Dataset name for structural connectome (e.g., 'Cook2019')
        synapse_type: 'chemical', 'electrical', or 'both'
        mode: 'binary' or 'weighted'
        symmetrize_electrical: Whether to symmetrize electrical synapse data

    Returns:
        dict with keys:
            'constrained': The constrained molecular connectome
            'structural': The structural connectome used
            'metadata': Dict with constraint parameters

    Example:
        >>> molecular = pu.assemble_nt_connectome('dopamine', ...)['binary']
        >>> result = constrain_molecular_connectome(molecular, 'Cook2019', 'chemical')
        >>> constrained = result['constrained']
    """
    neuron_order = molecular.index.tolist()

    # Load structural connectome
    structural = load_structural_connectome(
        synapse_type=synapse_type,
        dataset=structural_dataset,
        neuron_order=neuron_order,
        symmetrize_electrical=symmetrize_electrical
    )

    # Apply constraint
    constrained = apply_structural_constraint(molecular, structural, mode)

    # Build metadata
    metadata = {
        'structural_dataset': structural_dataset,
        'synapse_type': synapse_type,
        'mode': mode,
        'symmetrize_electrical': symmetrize_electrical
    }

    return {
        'constrained': constrained,
        'structural': structural,
        'metadata': metadata
    }


def get_available_structural_datasets() -> dict:
    """
    List available structural connectome datasets.

    Returns:
        dict with keys 'chemical_synapse' and 'electrical_synapse',
        each containing a list of available dataset names.

    Example:
        >>> datasets = get_available_structural_datasets()
        >>> print(datasets['chemical_synapse'])
        ['Cook2019', 'Varshney2011', ...]
    """
    assets = _load_assets()
    structural = assets['structural_connectomes']['preassembled']

    return {
        'chemical_synapse': list(structural['chemical_synapse'].keys()),
        'electrical_synapse': list(structural['electrical_synapse'].keys())
    }


# =============================================================================
# Comparison utilities
# =============================================================================

def compare_constrained_vs_unconstrained(
    molecular: pd.DataFrame,
    constrained: pd.DataFrame,
    threshold: float = 0
) -> dict:
    """
    Compare molecular connectome before and after structural constraint.

    Args:
        molecular: Original molecular connectome
        constrained: Structurally constrained molecular connectome
        threshold: Connection threshold (default 0, meaning > 0)

    Returns:
        dict with comparison statistics:
            'original_connections': Number of connections in original
            'constrained_connections': Number after constraint
            'removed_connections': Number removed by constraint
            'retention_rate': Fraction of connections retained
            'removed_pairs': List of (source, target) pairs that were removed
    """
    mol_mask = (molecular > threshold) & molecular.notna()
    con_mask = (constrained > threshold) & constrained.notna()

    n_original = mol_mask.sum().sum()
    n_constrained = con_mask.sum().sum()
    n_removed = n_original - n_constrained

    # Find removed pairs
    removed_mask = mol_mask & ~con_mask
    removed_pairs = []
    for source in molecular.index:
        for target in molecular.columns:
            if removed_mask.loc[source, target]:
                removed_pairs.append((source, target))

    return {
        'original_connections': int(n_original),
        'constrained_connections': int(n_constrained),
        'removed_connections': int(n_removed),
        'retention_rate': round(n_constrained / n_original, 4) if n_original > 0 else 0.0,
        'removed_pairs': removed_pairs
    }
