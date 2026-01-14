import pandas as pd


def get_release_vector(
    neurotransmitter: str,           # 'acetylcholine', 'dopamine', etc.
    markers: list[str],              # ['release'], ['synthesis', 'release'], ['eat-4'], etc.
    sources: list[str] = None,       # ['literature', 'reporter'] or None for all
    neuron_order: list = AllHermNeurons
) -> pd.Series:
    """
    Returns neuron vector. Value = 1 if ALL markers positive (AND gate), NaN if missing.
    Markers can be functional categories (mapped via gene_info) or specific gene names.
    """
    return None

def get_receptor_matrix(
    neurotransmitter: str,           # Filter receptors by their ligand
    sources: list[str],              # ['sequencing'], ['sequencing', 'reporter'], etc.
    gate: str = 'or',                # 'and' or 'or' across sources
    receptor_type: str = 'all',      # 'all', 'ionotropic', 'metabotropic'
    neuron_order: list = AllHermNeurons
) -> pd.DataFrame:
    """
    Returns neuron × receptor DataFrame. Values: 1/0/NaN.
    """
    return None

def assemble_nt_connectome(
    neurotransmitter: str,
    release_markers: list[str],
    release_sources: list[str] = None,
    receptor_sources: list[str] = ['sequencing'],
    receptor_gate: str = 'or',
    receptor_type: str = 'all',
    output_format: str = 'binary',   # 'per_pair', 'count', 'binary'
    neuron_order: list = AllHermNeurons
) -> pd.DataFrame | dict[str, pd.DataFrame]:
    """
    Returns:
    - 'per_pair': dict of {receptor_name: neuron × neuron DataFrame}
    - 'count': neuron × neuron DataFrame with receptor counts
    - 'binary': neuron × neuron DataFrame with 1 if count >= 1, else 0
    """
    return None