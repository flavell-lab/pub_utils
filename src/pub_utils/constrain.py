"""
DRAFT - Structural Constraints on Molecular Connectomes
========================================================

This module applies structural constraints to molecular connectomes.
A neuron can only receive a signal if it's structurally connected to the source.

Core logic:
    constrained[source, target] = molecular[source, target] AND structural[source, target]

Key design decisions:
    - Only chemical synapses are used to constrain neurotransmitter connectomes
      (directional signaling; electrical synapses are bidirectional gap junctions)
    - Default structural constraint = None (no constraint unless specified)
    - Modes: 'binary' (mask), 'weighted' (multiply)
"""

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from .constants import AllHermNeurons

# Default base path points to repository root; can be overridden
_DEFAULT_BASE_PATH = Path(__file__).parent.parent.parent


def _load_assets(base_path: Path = None) -> dict:
    """Load assets.json from base_path."""
    if base_path is None:
        base_path = _DEFAULT_BASE_PATH
    with open(base_path / "assets.json") as f:
        return json.load(f)


def _get_path(asset_path: str, base_path: Path = None) -> Path:
    """Convert relative asset path to absolute path."""
    if base_path is None:
        base_path = _DEFAULT_BASE_PATH
    return base_path / asset_path


def load_structural_connectome(
    dataset: str,
    neuron_order: list = None,
    base_path: Path = None,
) -> pd.DataFrame:
    """
    Load a chemical synapse structural connectome from preassembled files.

    Note: Only chemical synapses are used for constraining neurotransmitter
    connectomes (directional signaling).

    Args:
        dataset: Dataset name (e.g., 'Cook2019', 'Varshney2011', 'WhiteWhole')
        neuron_order: Canonical neuron ordering for output. Default: AllHermNeurons
        base_path: Repository root path. Default: auto-detected from module location.

    Returns:
        pd.DataFrame (source x target) with chemical synapse counts.

    Example:
        >>> structural = load_structural_connectome('Cook2019')
    """
    if neuron_order is None:
        neuron_order = AllHermNeurons

    assets = _load_assets(base_path)

    chem_assets = assets['structural_connectomes']['preassembled']['chemical_synapse']
    chem_path = chem_assets.get(dataset)
    if chem_path is None:
        available = list(chem_assets.keys())
        raise ValueError(
            f"No chemical synapse data for dataset '{dataset}'. "
            f"Available: {available}"
        )
    chem = pd.read_csv(_get_path(chem_path, base_path), index_col=0)
    chem = chem.reindex(index=neuron_order, columns=neuron_order).fillna(0)

    return chem


def get_available_structural_datasets(base_path: Path = None) -> list:
    """
    List available chemical synapse structural connectome datasets.

    Args:
        base_path: Repository root path. Default: auto-detected.

    Returns:
        List of available dataset names for chemical synapse connectomes.
    """
    assets = _load_assets(base_path)
    structural = assets['structural_connectomes']['preassembled']
    return list(structural['chemical_synapse'].keys())


def _compute_percentile_weights(structural: pd.DataFrame) -> pd.DataFrame:
    """
    Convert synapse counts to 0-5 graded weights based on percentiles.

    0 = no synapses
    1 = bottom 20% of non-zero values
    2 = 20-40 percentile
    3 = 40-60 percentile
    4 = 60-80 percentile
    5 = top 20 percentile
    """
    # Get all non-zero synapse counts
    nonzero_values = structural.values[structural.values > 0]

    if len(nonzero_values) == 0:
        return pd.DataFrame(0.0, index=structural.index, columns=structural.columns)

    # Compute percentile thresholds
    percentiles = [20, 40, 60, 80]
    thresholds = [np.percentile(nonzero_values, p) for p in percentiles]

    # Assign grades: 0 for zero, 1-5 based on percentile bin
    weights = pd.DataFrame(0.0, index=structural.index, columns=structural.columns)

    for i, j in zip(*np.where(structural.values > 0)):
        val = structural.values[i, j]
        if val <= thresholds[0]:
            grade = 1
        elif val <= thresholds[1]:
            grade = 2
        elif val <= thresholds[2]:
            grade = 3
        elif val <= thresholds[3]:
            grade = 4
        else:
            grade = 5
        weights.iloc[i, j] = grade

    return weights


def apply_structural_constraint(
    molecular: pd.DataFrame,
    structural: pd.DataFrame,
    mode: str = 'binary'
) -> pd.DataFrame:
    """
    Constrain a single molecular connectome by structural connectivity.

    Args:
        molecular: Molecular connectome (source x target)
        structural: Structural connectome (source x target)
        mode:
            'binary' - Mask: 1 if structural connection exists, else 0
            'weighted' - Graded 0-5 based on synapse count percentiles:
                         0=none, 1=bottom 20%, 2=20-40%, 3=40-60%, 4=60-80%, 5=top 20%

    Returns:
        Constrained molecular connectome.
    """
    common_idx = molecular.index.intersection(structural.index)
    common_col = molecular.columns.intersection(structural.columns)

    if len(common_idx) < len(molecular.index):
        missing = set(molecular.index) - set(common_idx)
        warnings.warn(
            f"{len(missing)} neurons in molecular but not structural. "
            f"First 5: {list(missing)[:5]}. Treated as no structural connection."
        )

    mol_aligned = molecular.loc[common_idx, common_col]
    struct_aligned = structural.loc[common_idx, common_col]

    if mode == 'binary':
        mask = (struct_aligned > 0).astype(float)
        constrained = mol_aligned * mask
    elif mode == 'weighted':
        weights = _compute_percentile_weights(struct_aligned)
        constrained = mol_aligned * weights
    else:
        raise ValueError(f"Unknown mode '{mode}'. Expected 'binary' or 'weighted'.")

    result = pd.DataFrame(0.0, index=molecular.index, columns=molecular.columns)
    result.loc[common_idx, common_col] = constrained
    return result


def constrain_assembly(
    assembly: dict,
    structural_dataset: str = 'Cook2019',
    mode: str = 'binary',
    base_path: Path = None,
) -> dict:
    """
    Apply structural constraint to a full assembly result.

    Constrains 'binary', 'count', and all matrices in 'per_pair'.

    Args:
        assembly: Output from assemble_nt_connectome() or assemble_npp_connectome()
                  with keys 'binary', 'count', 'per_pair'
        structural_dataset: Dataset name (e.g., 'Cook2019')
        mode: 'binary' or 'weighted'
        base_path: Repository root path. Default: auto-detected.

    Returns:
        dict with same structure as input, all matrices constrained:
            'binary': constrained binary matrix
            'count': constrained count matrix
            'per_pair': dict of constrained per-receptor matrices
            'structural': the structural connectome used
            'metadata': constraint parameters

    Example:
        >>> assembly = pu.assemble_nt_connectome('dopamine', ...)
        >>> constrained = constrain_assembly(assembly, 'Cook2019')
    """
    neuron_order = assembly['binary'].index.tolist()

    structural = load_structural_connectome(
        dataset=structural_dataset,
        neuron_order=neuron_order,
        base_path=base_path,
    )

    result = {
        'binary': apply_structural_constraint(assembly['binary'], structural, mode),
        'count': apply_structural_constraint(assembly['count'], structural, mode),
        'per_pair': {},
        'structural': structural,
        'metadata': {
            'structural_dataset': structural_dataset,
            'synapse_type': 'chemical',
            'mode': mode,
        }
    }

    for receptor_name, matrix in assembly['per_pair'].items():
        result['per_pair'][receptor_name] = apply_structural_constraint(
            matrix, structural, mode
        )

    return result


# For comparing constrained vs unconstrained connectomes, use:
#   from pub_utils.io import compare_connectomes
# Save both to CSV, then:
#   compare_connectomes('unconstrained.csv', 'constrained.csv')
