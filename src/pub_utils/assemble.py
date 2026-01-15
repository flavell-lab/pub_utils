import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from .constants import AllHermNeurons

# =============================================================================
# Module-level caches
# =============================================================================
_ASSETS_CACHE = None
_GENE_INFO_CACHE = None
_PAIRING_INFO_CACHE = {}

# Base path points to repository root (src/pub_utils/assemble.py -> repo root)
_BASE_PATH = Path(__file__).parent.parent.parent


# =============================================================================
# Data Loading Infrastructure (Phase 1)
# =============================================================================

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


def _load_gene_info() -> pd.DataFrame:
    """
    Load NT_uptake_synthesis_release_gene_info.csv.

    Returns DataFrame with columns: Neurotransmitter, UptakeGene, SynthesisGene, ReleaseGene
    Neurotransmitter names are lowercased for consistent lookup.
    """
    global _GENE_INFO_CACHE
    if _GENE_INFO_CACHE is None:
        assets = _load_assets()
        path = _get_path(assets["gene_info"]["neurotransmitter"])
        df = pd.read_csv(path)
        df["Neurotransmitter"] = df["Neurotransmitter"].str.lower()
        _GENE_INFO_CACHE = df
    return _GENE_INFO_CACHE


def _load_pairing_info(molecule_type: str = "neurotransmitter", source: str = None) -> pd.DataFrame:
    """
    Load receptor-ligand pairing info.

    Args:
        molecule_type: 'neurotransmitter' or 'neuropeptide'
        source: For NPP, one of 'Altun2013', 'Bentley2016', 'RipollSanchez2023'
                Ignored for neurotransmitter (single source).

    Returns:
        DataFrame with columns including: receptor, ligand, confidence, channel, GPCR
    """
    cache_key = (molecule_type, source)
    if cache_key not in _PAIRING_INFO_CACHE:
        assets = _load_assets()
        if molecule_type == "neurotransmitter":
            path = _get_path(assets["pairing_info"]["neurotransmitter"])
        else:
            if source is None:
                raise ValueError("source required for neuropeptide pairing info")
            path = _get_path(assets["pairing_info"]["neuropeptide"][source])
        df = pd.read_csv(path)
        df["ligand"] = df["ligand"].str.lower()
        _PAIRING_INFO_CACHE[cache_key] = df
    return _PAIRING_INFO_CACHE[cache_key]


def _load_release_data(
    molecule_type: str,
    source_method: str,
    source_dataset: str
) -> pd.DataFrame:
    """
    Load release expression matrix (neuronID × gene/neuropeptide).

    Args:
        molecule_type: 'neurotransmitter' or 'neuropeptide'
        source_method: 'literature', 'reporter', 'staining', 'sequencing'
        source_dataset: 'Bentley2016', 'Wang2024', 'RipollSanchez2023', etc.

    Returns:
        DataFrame with neuronID as index, genes/NPPs as columns, values 0/1.
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


def _load_receptor_data(
    molecule_type: str,
    neurotransmitter: str,
    receptor_type: str,
    source_method: str,
    source_dataset: str
) -> pd.DataFrame:
    """
    Load receptor expression matrix (neuronID × receptor).

    Args:
        molecule_type: 'neurotransmitter' or 'neuropeptide'
        neurotransmitter: NT name for filtering (e.g., 'dopamine'), or 'all'
        receptor_type: 'ionotropic', 'metabotropic', 'all'
        source_method: 'sequencing', 'reporter', 'literature'
        source_dataset: 'Fenyves2020', 'HobertLab', 'Bentley2016', etc.

    Returns:
        DataFrame with neuronID as index, receptor names as columns, values 0/1.
    """
    assets = _load_assets()
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


def _get_available_release_sources(molecule_type: str) -> list[str]:
    """
    Get all available release data sources for a molecule type.

    Returns list of source keys in format 'method:dataset'.
    """
    assets = _load_assets()
    sources = []

    release_assets = assets["release"].get(molecule_type, {})
    for method, datasets in release_assets.items():
        for dataset in datasets:
            sources.append(f"{method}:{dataset}")

    return sources


def _resolve_markers(neurotransmitter: str, markers: list[str]) -> list[str]:
    """
    Resolve functional category names to gene names.

    Args:
        neurotransmitter: NT name (e.g., 'dopamine')
        markers: List of functional categories or gene names
                 Categories: 'uptake', 'synthesis', 'release'
                 Gene names pass through unchanged.

    Returns:
        List of resolved gene names.

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


# =============================================================================
# Public API Functions
# =============================================================================

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
                 If None, uses all available sources.
        neuron_order: Canonical neuron ordering for output.

    Returns:
        pd.Series indexed by neuron_order.
        Values: 1 if ALL markers positive (AND gate), 0 if any marker 0, NaN if missing.

    Logic:
        1. Resolve functional marker names to gene names
        2. Load release matrices from specified sources
        3. For each source, extract columns matching resolved genes
        4. Apply AND gate across genes within each source
        5. Apply OR gate across sources (any positive source = positive)
        6. Align to neuron_order, fill missing with NaN
    """
    # Step 1: Resolve markers to gene names
    genes = _resolve_markers(neurotransmitter, markers)

    # Step 2: Determine sources
    if sources is None:
        sources = _get_available_release_sources("neurotransmitter")

    # Step 3: Load and combine data from each source
    source_results = []

    for source_key in sources:
        parts = source_key.split(":")
        if len(parts) != 2:
            raise ValueError(f"Invalid source format '{source_key}', expected 'method:dataset'")
        method, dataset = parts

        # Check if source exists in assets
        assets = _load_assets()
        release_assets = assets.get("release", {}).get("neurotransmitter", {})
        if method not in release_assets or dataset not in release_assets.get(method, {}):
            warnings.warn(f"Release source '{source_key}' not found in assets, skipping")
            continue

        df = _load_release_data("neurotransmitter", method, dataset)

        # Check which genes are available in this dataset
        available_genes = [g for g in genes if g in df.columns]
        if not available_genes:
            warnings.warn(
                f"Source '{source_key}' has no columns matching genes {genes}, skipping"
            )
            continue

        # AND gate across genes for this source (min of 0/1 values)
        # If any gene is 0, result is 0. If all genes are 1, result is 1.
        source_result = df[available_genes].min(axis=1)
        source_results.append(source_result)

    # Step 4: OR gate across sources
    if not source_results:
        # No data found, return NaN series
        return pd.Series(np.nan, index=neuron_order)

    # Combine all sources: OR gate = max across sources
    combined = pd.concat(source_results, axis=1)
    result = combined.max(axis=1)

    # Step 5: Align to neuron_order
    result = result.reindex(neuron_order)

    return result

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
        sources: List of source keys ['sequencing:Fenyves2020', 'literature:Bentley2016']
        gate: 'and' (require all sources) or 'or' (any source sufficient)
        receptor_type: 'all', 'ionotropic', or 'metabotropic'
        neuron_order: Canonical neuron ordering for rows.

    Returns:
        pd.DataFrame with neuron_order as index, receptor names as columns.
        Values: 1/0/NaN
    """
    # Step 1: Get valid receptors for this NT from pairing info
    pairing = _load_pairing_info("neurotransmitter")
    nt_pairs = pairing[pairing["ligand"] == neurotransmitter.lower()]

    # Filter by confidence >= 1
    nt_pairs = nt_pairs[nt_pairs["confidence"] >= 1]

    # Filter by receptor type
    if receptor_type == "ionotropic":
        nt_pairs = nt_pairs[nt_pairs["channel"] == 1]
    elif receptor_type == "metabotropic":
        nt_pairs = nt_pairs[nt_pairs["GPCR"] == 1]

    valid_receptors = set(nt_pairs["receptor"].unique())

    if not valid_receptors:
        warnings.warn(f"No valid receptors found for {neurotransmitter} with type={receptor_type}")
        return pd.DataFrame(index=neuron_order)

    # Step 2: Load receptor data from each source
    receptor_data_by_source = {}
    assets = _load_assets()

    for source_key in sources:
        parts = source_key.split(":")
        if len(parts) != 2:
            raise ValueError(f"Invalid source format '{source_key}', expected 'method:dataset'")
        method, dataset = parts

        # Try loading with specific receptor_type, fall back to 'all'
        df = None
        for rtype in [receptor_type, 'all']:
            # Check if this path exists in assets
            receptor_assets = assets.get("receptor", {}).get("neurotransmitter", {})

            # Check NT-specific first, then 'all' NT
            for nt_key in [neurotransmitter.lower(), 'all']:
                if nt_key not in receptor_assets:
                    continue
                nt_assets = receptor_assets[nt_key]
                if rtype not in nt_assets:
                    continue
                type_assets = nt_assets[rtype]
                if method not in type_assets:
                    continue
                method_assets = type_assets[method]
                if dataset not in method_assets:
                    continue

                # Found a valid path, load it
                df = _load_receptor_data("neurotransmitter", nt_key, rtype, method, dataset)
                break

            if df is not None:
                break

        if df is None:
            warnings.warn(f"Receptor source '{source_key}' not found for {neurotransmitter}, skipping")
            continue

        # Filter to valid receptors that exist in this dataset
        available = [r for r in valid_receptors if r in df.columns]
        if available:
            receptor_data_by_source[source_key] = df[available]
        else:
            warnings.warn(
                f"Source '{source_key}' has no columns matching valid receptors for {neurotransmitter}"
            )

    if not receptor_data_by_source:
        warnings.warn(f"No receptor data loaded for {neurotransmitter}")
        return pd.DataFrame(index=neuron_order)

    # Step 3: Combine across sources per receptor
    all_receptors = set()
    for df in receptor_data_by_source.values():
        all_receptors.update(df.columns)

    result = pd.DataFrame(index=neuron_order, columns=sorted(all_receptors), dtype=float)

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

def assemble_nt_connectome(
    neurotransmitter: str,
    release_markers: list[str],
    release_sources: list[str] = None,
    receptor_sources: list[str] = None,
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
        release_sources: Sources for release data (e.g., ['literature:Bentley2016']).
                        None = use all available.
        receptor_sources: Sources for receptor data (e.g., ['reporter:Muralidhara2025']).
                         Required - no default to force explicit choice.
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
    if receptor_sources is None:
        raise ValueError("receptor_sources must be specified explicitly")

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
                np.nan, index=neuron_order, columns=neuron_order, dtype=float
            )

    # Step 3: Compute per-receptor connectomes via outer product
    per_pair = {}

    for receptor_name in receptor.columns:
        receptor_vec = receptor[receptor_name]

        # Outer product: release[source] × receptor[target]
        # Connection exists if source releases AND target has receptor
        release_arr = release.values.reshape(-1, 1)
        receptor_arr = receptor_vec.values.reshape(1, -1)

        # Multiply handles NaN propagation: NaN * x = NaN
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
        # First get count matrix, then threshold
        count_matrix = assemble_nt_connectome(
            neurotransmitter, release_markers, release_sources,
            receptor_sources, receptor_gate, receptor_type,
            'count', neuron_order
        )
        binary_matrix = (count_matrix >= 1).astype(float)
        binary_matrix[count_matrix.isna()] = np.nan
        return binary_matrix


# =============================================================================
# Neuropeptide (NPP) Assembly Functions
# =============================================================================

def _load_npp_pairing_info(source: str) -> pd.DataFrame:
    """
    Load NPP receptor-ligand pairing info.

    All sources now have standardized columns: 'receptor', 'ligand'

    Returns DataFrame with lowercase receptor/ligand columns.
    """
    assets = _load_assets()
    path = _get_path(assets["pairing_info"]["neuropeptide"][source])
    pairing = pd.read_csv(path)

    # Drop rows with empty receptor/ligand
    pairing = pairing.dropna(subset=["receptor", "ligand"])

    # Lowercase for consistent matching
    pairing["receptor"] = pairing["receptor"].str.lower()
    pairing["ligand"] = pairing["ligand"].str.lower()

    return pairing


def _load_npp_release_data(source_method: str, source_dataset: str) -> pd.DataFrame:
    """Load NPP release expression matrix (neuronID × neuropeptide)."""
    return _load_release_data("neuropeptide", source_method, source_dataset)


def _load_npp_receptor_data(source_method: str, source_dataset: str) -> pd.DataFrame:
    """
    Load NPP receptor expression matrix (neuronID × receptor).

    Navigates assets.json structure: receptor.neuropeptide.metabotropic.{method}.{dataset}
    """
    assets = _load_assets()

    # NPP receptors are under receptor.neuropeptide.metabotropic
    receptor_assets = assets.get("receptor", {}).get("neuropeptide", {}).get("metabotropic", {})

    if source_method not in receptor_assets:
        raise ValueError(f"No NPP receptor data for method '{source_method}'")
    if source_dataset not in receptor_assets[source_method]:
        raise ValueError(f"No NPP receptor data for {source_method}/{source_dataset}")

    path_str = receptor_assets[source_method][source_dataset]
    df = pd.read_csv(_get_path(path_str))
    df = df.set_index("neuronID")

    # Lowercase column names for consistent matching
    df.columns = df.columns.str.lower()

    return df


def get_npp_release_vector(
    neuropeptide: str,
    sources: list[str] = None,
    neuron_order: list = AllHermNeurons
) -> pd.Series:
    """
    Returns neuron vector with release status for a neuropeptide.

    Args:
        neuropeptide: NPP name (e.g., 'flp-1', 'nlp-12')
        sources: List of source keys ['literature:Bentley2016', 'sequencing:RipollSanchez2023']
                 If None, uses all available sources.
        neuron_order: Canonical neuron ordering for output.

    Returns:
        pd.Series indexed by neuron_order. Values: 1/0/NaN
    """
    npp_lower = neuropeptide.lower()

    # Determine sources
    if sources is None:
        sources = _get_available_release_sources("neuropeptide")

    source_results = []
    assets = _load_assets()

    for source_key in sources:
        parts = source_key.split(":")
        if len(parts) != 2:
            raise ValueError(f"Invalid source format '{source_key}', expected 'method:dataset'")
        method, dataset = parts

        # Check if source exists
        release_assets = assets.get("release", {}).get("neuropeptide", {})
        if method not in release_assets or dataset not in release_assets.get(method, {}):
            warnings.warn(f"NPP release source '{source_key}' not found, skipping")
            continue

        df = _load_npp_release_data(method, dataset)
        df.columns = df.columns.str.lower()

        if npp_lower not in df.columns:
            continue  # This NPP not in this dataset

        source_results.append(df[npp_lower])

    if not source_results:
        warnings.warn(f"No release data found for neuropeptide '{neuropeptide}'")
        return pd.Series(np.nan, index=neuron_order)

    # OR gate across sources
    combined = pd.concat(source_results, axis=1)
    result = combined.max(axis=1)

    return result.reindex(neuron_order)


def get_npp_receptor_matrix(
    neuropeptide: str,
    sources: list[str],
    pairing_source: str = "RipollSanchez2023",
    gate: str = 'or',
    neuron_order: list = AllHermNeurons
) -> pd.DataFrame:
    """
    Returns neuron × receptor matrix for receptors of specified neuropeptide.

    Args:
        neuropeptide: NPP name to filter receptors by ligand
        sources: List of source keys ['sequencing:RipollSanchez2023', 'literature:Bentley2016']
        pairing_source: Which pairing info to use ('Altun2013', 'Bentley2016', 'RipollSanchez2023')
        gate: 'and' or 'or' across sources
        neuron_order: Canonical neuron ordering for rows.

    Returns:
        pd.DataFrame with neuron_order as index, receptor names as columns.
    """
    npp_lower = neuropeptide.lower()

    # Get valid receptors for this NPP from pairing info
    pairing = _load_npp_pairing_info(pairing_source)
    npp_pairs = pairing[pairing["ligand"] == npp_lower]

    valid_receptors = set(npp_pairs["receptor"].unique())

    if not valid_receptors:
        warnings.warn(f"No valid receptors found for {neuropeptide} in {pairing_source} pairing info")
        return pd.DataFrame(index=neuron_order)

    # Load receptor data from each source
    receptor_data_by_source = {}

    for source_key in sources:
        parts = source_key.split(":")
        if len(parts) != 2:
            raise ValueError(f"Invalid source format '{source_key}', expected 'method:dataset'")
        method, dataset = parts

        try:
            df = _load_npp_receptor_data(method, dataset)
        except ValueError as e:
            warnings.warn(f"NPP receptor source '{source_key}' not found: {e}")
            continue

        # Filter to valid receptors that exist in this dataset
        available = [r for r in valid_receptors if r in df.columns]
        if available:
            receptor_data_by_source[source_key] = df[available]

    if not receptor_data_by_source:
        warnings.warn(f"No receptor data loaded for {neuropeptide}")
        return pd.DataFrame(index=neuron_order)

    # Combine across sources per receptor
    all_receptors = set()
    for df in receptor_data_by_source.values():
        all_receptors.update(df.columns)

    result = pd.DataFrame(index=neuron_order, columns=sorted(all_receptors), dtype=float)

    for receptor in all_receptors:
        receptor_values = []
        for source_key, df in receptor_data_by_source.items():
            if receptor in df.columns:
                receptor_values.append(df[receptor])

        if not receptor_values:
            continue

        combined = pd.concat(receptor_values, axis=1)

        if gate == 'or':
            result[receptor] = combined.max(axis=1).reindex(neuron_order)
        else:
            result[receptor] = combined.min(axis=1).reindex(neuron_order)

    return result


def assemble_npp_connectome(
    neuropeptide: str = None,
    release_sources: list[str] = None,
    receptor_sources: list[str] = None,
    pairing_source: str = "RipollSanchez2023",
    receptor_gate: str = 'or',
    output_format: str = 'binary',
    neuron_order: list = AllHermNeurons
) -> pd.DataFrame | dict[str, pd.DataFrame]:
    """
    Assemble a molecular connectome for a neuropeptide.

    Args:
        neuropeptide: NPP name (e.g., 'flp-1'). Required.
        release_sources: Sources for release data. None = use all available.
        receptor_sources: Sources for receptor data. Required.
        pairing_source: Which pairing info ('Altun2013', 'Bentley2016', 'RipollSanchez2023')
        receptor_gate: 'or' or 'and' across receptor sources
        output_format:
            'per_pair' - dict of {receptor: source×target DataFrame}
            'count' - source×target DataFrame with receptor counts
            'binary' - source×target DataFrame with 1 if any connection
        neuron_order: Canonical neuron ordering

    Returns:
        DataFrame or dict depending on output_format.
        Matrix semantics: rows = source neurons, columns = target neurons.
    """
    if neuropeptide is None:
        raise ValueError("neuropeptide must be specified")
    if receptor_sources is None:
        raise ValueError("receptor_sources must be specified explicitly")

    # Step 1: Get release vector (source capability)
    release = get_npp_release_vector(neuropeptide, release_sources, neuron_order)

    # Step 2: Get receptor matrix (target capability)
    receptor = get_npp_receptor_matrix(
        neuropeptide, receptor_sources, pairing_source, receptor_gate, neuron_order
    )

    if receptor.empty:
        if output_format == 'per_pair':
            return {}
        else:
            return pd.DataFrame(
                np.nan, index=neuron_order, columns=neuron_order, dtype=float
            )

    # Step 3: Compute per-receptor connectomes via outer product
    per_pair = {}

    for receptor_name in receptor.columns:
        receptor_vec = receptor[receptor_name]

        release_arr = release.values.reshape(-1, 1)
        receptor_arr = receptor_vec.values.reshape(1, -1)

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
        count_matrix = pd.DataFrame(
            0.0, index=neuron_order, columns=neuron_order
        )
        for receptor_name, matrix in per_pair.items():
            count_matrix = count_matrix.add(matrix.fillna(0))

        all_nan_mask = pd.DataFrame(True, index=neuron_order, columns=neuron_order)
        for matrix in per_pair.values():
            all_nan_mask = all_nan_mask & matrix.isna()
        count_matrix[all_nan_mask] = np.nan

        return count_matrix

    else:  # 'binary'
        count_matrix = assemble_npp_connectome(
            neuropeptide, release_sources, receptor_sources,
            pairing_source, receptor_gate, 'count', neuron_order
        )
        binary_matrix = (count_matrix >= 1).astype(float)
        binary_matrix[count_matrix.isna()] = np.nan
        return binary_matrix


# =============================================================================
# Connectome I/O with Metadata
# =============================================================================

def save_connectome(
    connectome: pd.DataFrame,
    filepath: str,
    metadata: dict
) -> None:
    """
    Save a connectome CSV with a JSON sidecar containing metadata.

    Args:
        connectome: The connectome DataFrame to save
        filepath: Path to save CSV (e.g., 'connectomes/dk_assembly/dopamine_01.csv')
        metadata: Dict with assembly parameters. Recommended keys:
            - molecule: str (e.g., 'dopamine', 'flp-1')
            - molecule_type: str ('neurotransmitter' or 'neuropeptide')
            - release_markers: list[str] (for NT only)
            - release_sources: list[str]
            - receptor_sources: list[str]
            - receptor_gate: str
            - receptor_type: str (for NT only)
            - pairing_source: str (for NPP only)
            - output_format: str
            - description: str (optional)

    Creates:
        - {filepath} - CSV file with connectome data
        - {filepath}.json - JSON file with metadata
    """
    from datetime import datetime

    # Save CSV
    connectome.to_csv(filepath)

    # Add automatic metadata
    full_metadata = {
        "created": datetime.now().isoformat(),
        "shape": list(connectome.shape),
        "total_connections": int(connectome.sum().sum()),
        "nonzero_entries": int((connectome > 0).sum().sum()),
        **metadata
    }

    # Save JSON sidecar
    json_path = filepath + ".json"
    with open(json_path, "w") as f:
        json.dump(full_metadata, f, indent=2)


def load_connectome(filepath: str) -> tuple[pd.DataFrame, dict]:
    """
    Load a connectome CSV and its JSON metadata sidecar.

    Args:
        filepath: Path to CSV file

    Returns:
        Tuple of (connectome DataFrame, metadata dict)
        If no JSON sidecar exists, metadata will be empty dict.
    """
    connectome = pd.read_csv(filepath, index_col=0)

    json_path = filepath + ".json"
    metadata = {}
    if Path(json_path).exists():
        with open(json_path) as f:
            metadata = json.load(f)

    return connectome, metadata


# =============================================================================
# Structural Constraints
# =============================================================================

def load_structural_connectome(
    synapse_type: str,
    dataset: str,
    neuron_order: list = AllHermNeurons
) -> pd.DataFrame:
    """
    Load a structural connectome from preassembled files.

    Args:
        synapse_type: 'chemical', 'electrical', or 'both'
        dataset: Dataset name (e.g., 'Cook2019', 'Varshney2011', 'WhiteWhole')
        neuron_order: Canonical neuron ordering for output

    Returns:
        pd.DataFrame (source × target) with synapse counts.
        For 'both', returns union (sum) of chemical and electrical.
    """
    assets = _load_assets()

    if synapse_type in ['chemical', 'both']:
        chem_path = assets['structural_connectomes']['preassembled']['chemical_synapse'].get(dataset)
        if chem_path is None:
            raise ValueError(f"No chemical synapse data for dataset '{dataset}'")
        chem = pd.read_csv(_get_path(chem_path), index_col=0)
        chem = chem.reindex(index=neuron_order, columns=neuron_order).fillna(0)

    if synapse_type in ['electrical', 'both']:
        elec_path = assets['structural_connectomes']['preassembled']['electrical_synapse'].get(dataset)
        if elec_path is None:
            raise ValueError(f"No electrical synapse data for dataset '{dataset}'")
        elec = pd.read_csv(_get_path(elec_path), index_col=0)
        elec = elec.reindex(index=neuron_order, columns=neuron_order).fillna(0)

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
        molecular: Molecular connectome (source × target)
        structural: Structural connectome (source × target), e.g., from load_structural_connectome()
        mode:
            'binary' - Mask: keep molecular connection only if structural connection exists
            'weighted' - Multiply: weight molecular by structural synapse count

    Returns:
        Constrained molecular connectome with same shape as input.
    """
    # Align indices
    common_idx = molecular.index.intersection(structural.index)
    common_col = molecular.columns.intersection(structural.columns)

    if len(common_idx) < len(molecular.index):
        missing = set(molecular.index) - set(common_idx)
        warnings.warn(f"{len(missing)} neurons in molecular but not structural: {list(missing)[:5]}...")

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
        raise ValueError(f"Unknown mode '{mode}', expected 'binary' or 'weighted'")

    # Restore full shape with zeros for missing neurons
    result = pd.DataFrame(0.0, index=molecular.index, columns=molecular.columns)
    result.loc[common_idx, common_col] = constrained

    return result