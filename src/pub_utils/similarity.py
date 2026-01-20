"""
Neuron Similarity Module

Computes pairwise similarity between 302 C. elegans neurons across four dimensions:
1. NT/NPP Release - Which neurons release similar neurotransmitters/neuropeptides
2. NT/NPP Reception - Which neurons express similar receptors
3. Sensory Type - Which neurons respond to similar sensory modalities
4. Anatomical Process - Which neurons have processes in similar anatomical regions

All similarity matrices are 302x302 DataFrames with AllHermNeurons as index/columns.
Values range from 0 (no similarity) to 1 (identical profiles).
"""

import warnings
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .constants import AllHermNeurons
from .core import NeuronFeatures
from .assemble import (
    get_release_vector,
    get_receptor_matrix,
    get_npp_release_vector,
    get_npp_receptor_matrix,
    _load_assets,
    _get_path,
    _load_gene_info,
    _get_available_release_sources,
)

# Base path points to repository root
_BASE_PATH = Path(__file__).parent.parent.parent


# =============================================================================
# Profile Aggregation Functions
# =============================================================================

def get_sensory_profile(neuron_order: list = AllHermNeurons) -> pd.DataFrame:
    """
    Extract sensory type features from neuron_features.csv.

    Args:
        neuron_order: Canonical neuron ordering

    Returns:
        DataFrame (302 x 13 sensory types), values 0/1
        Columns are sensory modalities (chemical, odor, thermal, etc.)
    """
    assets = _load_assets()
    path = _get_path(assets["neuron_features"]["neuroanatomy"])
    df = pd.read_csv(path)

    nf = NeuronFeatures(df)

    # Get sensory type columns
    sensory_matrix = nf.get_category_features('sensoryType')
    sensory_names = nf.get_category_names('sensoryType')

    # Clean column names (remove 'sensoryType:' prefix)
    clean_names = [name.replace('sensoryType:', '') for name in sensory_names]

    result = pd.DataFrame(
        sensory_matrix,
        index=nf.neuron_ids,
        columns=clean_names
    )

    return result.reindex(neuron_order)


def get_process_profile(neuron_order: list = AllHermNeurons) -> pd.DataFrame:
    """
    Extract anatomical process features from neuron_features.csv.

    Args:
        neuron_order: Canonical neuron ordering

    Returns:
        DataFrame (302 x 28 process locations), values 0/1
        Columns are anatomical locations (nerve_ring, dorsal_cord, etc.)
    """
    assets = _load_assets()
    path = _get_path(assets["neuron_features"]["neuroanatomy"])
    df = pd.read_csv(path)

    nf = NeuronFeatures(df)

    # Get process columns
    process_matrix = nf.get_category_features('process')
    process_names = nf.get_category_names('process')

    # Clean column names (remove 'process:' prefix)
    clean_names = [name.replace('process:', '') for name in process_names]

    result = pd.DataFrame(
        process_matrix,
        index=nf.neuron_ids,
        columns=clean_names
    )

    return result.reindex(neuron_order)


def get_nt_release_profile(
    neurotransmitters: list[str] = None,
    sources: list[str] = None,
    gate: str = 'or',
    neuron_order: list = AllHermNeurons
) -> pd.DataFrame:
    """
    Build NT release profile matrix (neurons x neurotransmitters).

    Args:
        neurotransmitters: List of NTs to include. None = all available
                          ['acetylcholine', 'dopamine', 'gaba', 'glutamate',
                           'serotonin', 'tyramine', 'octopamine', 'betaine']
        sources: Release data sources (e.g., ['literature:Bentley2016'])
        gate: 'or' or 'and' across sources
        neuron_order: Canonical neuron ordering

    Returns:
        DataFrame (302 x n_NTs), values 0/1/NaN
    """
    if neurotransmitters is None:
        # Get all NTs from gene info
        gene_info = _load_gene_info()
        neurotransmitters = gene_info['Neurotransmitter'].unique().tolist()

    profiles = {}
    for nt in neurotransmitters:
        try:
            vec = get_release_vector(nt, ['release'], sources, gate, neuron_order)
            profiles[nt] = vec
        except (ValueError, KeyError) as e:
            warnings.warn(f"Could not get release vector for {nt}: {e}")
            profiles[nt] = pd.Series(np.nan, index=neuron_order)

    return pd.DataFrame(profiles, index=neuron_order)


def get_nt_receptor_profile(
    neurotransmitters: list[str] = None,
    sources: list[str] = None,
    gate: str = 'or',
    neuron_order: list = AllHermNeurons
) -> pd.DataFrame:
    """
    Build NT receptor profile matrix (neurons x receptors).

    Aggregates all receptors across specified neurotransmitters.

    Args:
        neurotransmitters: Filter receptors by ligand. None = all NTs
        sources: Receptor data sources
        gate: 'or' or 'and' across sources
        neuron_order: Canonical neuron ordering

    Returns:
        DataFrame (302 x n_receptors), values 0/1/NaN
    """
    if neurotransmitters is None:
        gene_info = _load_gene_info()
        neurotransmitters = gene_info['Neurotransmitter'].unique().tolist()

    if sources is None:
        # Default to Fenyves2020 sequencing + HobertLab reporter
        sources = ['sequencing:Fenyves2020', 'reporter:HobertLab']

    all_receptors = []
    for nt in neurotransmitters:
        try:
            receptor_df = get_receptor_matrix(nt, sources, gate, 'all', neuron_order)
            if not receptor_df.empty:
                all_receptors.append(receptor_df)
        except (ValueError, KeyError) as e:
            warnings.warn(f"Could not get receptor matrix for {nt}: {e}")

    if not all_receptors:
        return pd.DataFrame(index=neuron_order)

    # Concatenate all receptor matrices, handling overlapping columns
    result = pd.concat(all_receptors, axis=1)

    # Remove duplicate columns (same receptor from different NTs)
    result = result.loc[:, ~result.columns.duplicated()]

    return result


def get_npp_release_profile(
    neuropeptides: list[str] = None,
    sources: list[str] = None,
    gate: str = 'or',
    neuron_order: list = AllHermNeurons
) -> pd.DataFrame:
    """
    Build NPP release profile matrix (neurons x neuropeptides).

    Args:
        neuropeptides: List of NPPs to include. None = all available (~109)
        sources: Release data sources
        gate: 'or' or 'and' across sources
        neuron_order: Canonical neuron ordering

    Returns:
        DataFrame (302 x ~109 NPPs), values 0/1/NaN
    """
    if sources is None:
        sources = _get_available_release_sources('neuropeptide')

    if neuropeptides is None:
        # Load one source to get all available NPPs
        assets = _load_assets()
        # Use RipollSanchez2023 sequencing as it has the most NPPs
        path = _get_path(assets['release']['neuropeptide']['sequencing']['RipollSanchez2023'])
        npp_df = pd.read_csv(path)
        neuropeptides = [col for col in npp_df.columns if col != 'neuronID']

    profiles = {}
    for npp in neuropeptides:
        vec = get_npp_release_vector(npp, sources, gate, neuron_order)
        profiles[npp] = vec

    return pd.DataFrame(profiles, index=neuron_order)


def get_npp_receptor_profile(
    neuropeptides: list[str] = None,
    sources: list[str] = None,
    pairing_source: str = 'RipollSanchez2023',
    gate: str = 'or',
    neuron_order: list = AllHermNeurons
) -> pd.DataFrame:
    """
    Build NPP receptor profile matrix (neurons x receptors).

    Args:
        neuropeptides: Filter receptors by ligand. None = use all NPPs
        sources: Receptor data sources
        pairing_source: Which pairing info to use
        gate: 'or' or 'and' across sources
        neuron_order: Canonical neuron ordering

    Returns:
        DataFrame (302 x ~154 receptors), values 0/1/NaN
    """
    if sources is None:
        sources = ['sequencing:RipollSanchez2023']

    if neuropeptides is None:
        # Get all NPPs from release data
        assets = _load_assets()
        path = _get_path(assets['release']['neuropeptide']['sequencing']['RipollSanchez2023'])
        npp_df = pd.read_csv(path)
        neuropeptides = [col for col in npp_df.columns if col != 'neuronID']

    all_receptors = []
    for npp in neuropeptides:
        try:
            receptor_df = get_npp_receptor_matrix(npp, sources, pairing_source, gate, neuron_order)
            if not receptor_df.empty:
                all_receptors.append(receptor_df)
        except (ValueError, KeyError):
            pass  # Skip NPPs without receptor info

    if not all_receptors:
        return pd.DataFrame(index=neuron_order)

    # Concatenate and deduplicate
    result = pd.concat(all_receptors, axis=1)
    result = result.loc[:, ~result.columns.duplicated()]

    return result


# =============================================================================
# Core Similarity Computation
# =============================================================================

def _jaccard_similarity(a: np.ndarray, b: np.ndarray) -> float:
    """
    Compute Jaccard similarity between two binary vectors.

    Jaccard = |A intersection B| / |A union B|

    Handles NaN by excluding positions where either value is NaN.
    Returns NaN if no valid positions or union is empty.
    """
    # Mask NaN positions
    valid = ~(np.isnan(a) | np.isnan(b))

    if not np.any(valid):
        return np.nan

    a_valid = a[valid].astype(bool)
    b_valid = b[valid].astype(bool)

    intersection = np.sum(a_valid & b_valid)
    union = np.sum(a_valid | b_valid)

    if union == 0:
        # Both vectors are all zeros - define as similarity 1.0 (identical)
        return 1.0

    return intersection / union


def _cosine_similarity(a: np.ndarray, b: np.ndarray) -> float:
    """
    Compute cosine similarity between two vectors.

    Cosine = (A dot B) / (||A|| * ||B||)

    Handles NaN by excluding positions where either value is NaN.
    Returns NaN if no valid positions or either norm is zero.
    """
    # Mask NaN positions
    valid = ~(np.isnan(a) | np.isnan(b))

    if not np.any(valid):
        return np.nan

    a_valid = a[valid]
    b_valid = b[valid]

    dot_product = np.dot(a_valid, b_valid)
    norm_a = np.linalg.norm(a_valid)
    norm_b = np.linalg.norm(b_valid)

    if norm_a == 0 or norm_b == 0:
        # Zero vectors - if both zero, identical; otherwise, undefined
        if norm_a == 0 and norm_b == 0:
            return 1.0
        return 0.0

    return dot_product / (norm_a * norm_b)


def compute_pairwise_similarity(
    profile: pd.DataFrame,
    metric: str = 'jaccard',
    min_features: int = 1
) -> pd.DataFrame:
    """
    Compute pairwise similarity matrix from a feature profile.

    Args:
        profile: DataFrame (neurons x features), values 0/1/NaN
        metric: 'jaccard' (binary) or 'cosine' (continuous)
        min_features: Minimum valid features required for similarity
                     (pairs with fewer set to NaN)

    Returns:
        DataFrame (n_neurons x n_neurons) with similarity values [0, 1]
        Diagonal is 1.0 (self-similarity)
    """
    if metric == 'jaccard':
        sim_func = _jaccard_similarity
    elif metric == 'cosine':
        sim_func = _cosine_similarity
    else:
        raise ValueError(f"Unknown metric '{metric}', expected 'jaccard' or 'cosine'")

    neurons = profile.index.tolist()
    n_neurons = len(neurons)
    matrix = profile.values

    # Pre-allocate result matrix
    result = np.zeros((n_neurons, n_neurons))

    # Compute pairwise similarities
    for i in range(n_neurons):
        result[i, i] = 1.0  # Self-similarity
        for j in range(i + 1, n_neurons):
            sim = sim_func(matrix[i], matrix[j])
            result[i, j] = sim
            result[j, i] = sim  # Symmetric

    return pd.DataFrame(result, index=neurons, columns=neurons)


# =============================================================================
# High-Level API (Dimension-Specific)
# =============================================================================

def compute_release_similarity(
    molecule_type: str = 'all',
    metric: str = 'jaccard',
    sources: dict = None,
    gate: str = 'or',
    neuron_order: list = AllHermNeurons
) -> dict:
    """
    Compute neuron similarity based on neurotransmitter/neuropeptide release.

    Args:
        molecule_type: 'nt' (neurotransmitters only), 'npp' (neuropeptides only),
                      or 'all' (combined)
        metric: Similarity metric ('jaccard' or 'cosine')
        sources: Optional dict {'nt': [...], 'npp': [...]} of source keys
        gate: Gate for combining sources
        neuron_order: Neuron ordering

    Returns:
        dict with keys:
            'similarity': 302x302 DataFrame of similarity scores
            'profile': Feature profile used for computation
            'metadata': Parameters used
    """
    profiles = []

    if molecule_type in ['nt', 'all']:
        nt_sources = sources.get('nt') if sources else None
        nt_profile = get_nt_release_profile(
            neurotransmitters=None, sources=nt_sources, gate=gate, neuron_order=neuron_order
        )
        profiles.append(nt_profile)

    if molecule_type in ['npp', 'all']:
        npp_sources = sources.get('npp') if sources else None
        npp_profile = get_npp_release_profile(
            neuropeptides=None, sources=npp_sources, gate=gate, neuron_order=neuron_order
        )
        profiles.append(npp_profile)

    if not profiles:
        raise ValueError(f"Unknown molecule_type '{molecule_type}'")

    # Combine profiles
    combined_profile = pd.concat(profiles, axis=1)

    # Compute similarity
    similarity = compute_pairwise_similarity(combined_profile, metric=metric)

    return {
        'similarity': similarity,
        'profile': combined_profile,
        'metadata': {
            'molecule_type': molecule_type,
            'metric': metric,
            'gate': gate,
            'n_features': combined_profile.shape[1]
        }
    }


def compute_receptor_similarity(
    molecule_type: str = 'all',
    metric: str = 'jaccard',
    sources: dict = None,
    gate: str = 'or',
    neuron_order: list = AllHermNeurons
) -> dict:
    """
    Compute neuron similarity based on receptor expression.

    Args:
        molecule_type: 'nt', 'npp', or 'all'
        metric: Similarity metric ('jaccard' or 'cosine')
        sources: Optional dict {'nt': [...], 'npp': [...]} of source keys
        gate: Gate for combining sources
        neuron_order: Neuron ordering

    Returns:
        dict with 'similarity', 'profile', 'metadata'
    """
    profiles = []

    if molecule_type in ['nt', 'all']:
        nt_sources = sources.get('nt') if sources else None
        nt_profile = get_nt_receptor_profile(
            neurotransmitters=None, sources=nt_sources, gate=gate, neuron_order=neuron_order
        )
        if not nt_profile.empty:
            profiles.append(nt_profile)

    if molecule_type in ['npp', 'all']:
        npp_sources = sources.get('npp') if sources else None
        npp_profile = get_npp_receptor_profile(
            neuropeptides=None, sources=npp_sources, gate=gate, neuron_order=neuron_order
        )
        if not npp_profile.empty:
            profiles.append(npp_profile)

    if not profiles:
        raise ValueError(f"No receptor data found for molecule_type '{molecule_type}'")

    # Combine profiles
    combined_profile = pd.concat(profiles, axis=1)
    combined_profile = combined_profile.loc[:, ~combined_profile.columns.duplicated()]

    # Compute similarity
    similarity = compute_pairwise_similarity(combined_profile, metric=metric)

    return {
        'similarity': similarity,
        'profile': combined_profile,
        'metadata': {
            'molecule_type': molecule_type,
            'metric': metric,
            'gate': gate,
            'n_features': combined_profile.shape[1]
        }
    }


def compute_sensory_similarity(
    metric: str = 'jaccard',
    neuron_order: list = AllHermNeurons
) -> dict:
    """
    Compute neuron similarity based on sensory modality responses.

    Args:
        metric: Similarity metric ('jaccard' recommended for binary data)
        neuron_order: Neuron ordering

    Returns:
        dict with 'similarity', 'profile', 'metadata'
    """
    profile = get_sensory_profile(neuron_order)
    similarity = compute_pairwise_similarity(profile, metric=metric)

    return {
        'similarity': similarity,
        'profile': profile,
        'metadata': {
            'metric': metric,
            'n_features': profile.shape[1],
            'feature_names': profile.columns.tolist()
        }
    }


def compute_process_similarity(
    metric: str = 'jaccard',
    neuron_order: list = AllHermNeurons
) -> dict:
    """
    Compute neuron similarity based on anatomical process locations.

    Args:
        metric: Similarity metric ('jaccard' recommended for binary data)
        neuron_order: Neuron ordering

    Returns:
        dict with 'similarity', 'profile', 'metadata'
    """
    profile = get_process_profile(neuron_order)
    similarity = compute_pairwise_similarity(profile, metric=metric)

    return {
        'similarity': similarity,
        'profile': profile,
        'metadata': {
            'metric': metric,
            'n_features': profile.shape[1],
            'feature_names': profile.columns.tolist()
        }
    }


# =============================================================================
# Convenience Functions
# =============================================================================

def compute_all_similarities(
    neuron_order: list = AllHermNeurons,
    save_path: str = None
) -> dict:
    """
    Compute all four similarity dimensions at once.

    Args:
        neuron_order: Canonical neuron ordering
        save_path: Optional directory path to save CSV outputs

    Returns:
        dict with keys: 'release', 'receptor', 'sensory', 'process'
        Each value is output from the corresponding compute_*_similarity() function
    """
    results = {
        'sensory': compute_sensory_similarity(neuron_order=neuron_order),
        'process': compute_process_similarity(neuron_order=neuron_order),
        'release': compute_release_similarity(molecule_type='all', neuron_order=neuron_order),
        'receptor': compute_receptor_similarity(molecule_type='all', neuron_order=neuron_order),
    }

    if save_path:
        save_dir = Path(save_path)
        save_dir.mkdir(parents=True, exist_ok=True)

        for dim_name, result in results.items():
            filepath = save_dir / f"{dim_name}_similarity.csv"
            result['similarity'].to_csv(filepath)
            print(f"Saved {filepath}")

    return results


def find_similar_neurons(
    similarity_matrix: pd.DataFrame,
    neuron: str,
    top_n: int = 10,
    threshold: float = None
) -> pd.DataFrame:
    """
    Find neurons most similar to a given neuron.

    Args:
        similarity_matrix: Output from compute_*_similarity()['similarity']
        neuron: Query neuron ID
        top_n: Number of most similar neurons to return
        threshold: Optional minimum similarity threshold

    Returns:
        DataFrame with columns: neuron_id, similarity
        Sorted by similarity descending (excludes self)
    """
    if neuron not in similarity_matrix.index:
        raise ValueError(f"Neuron '{neuron}' not found in similarity matrix")

    similarities = similarity_matrix.loc[neuron].copy()

    # Remove self
    similarities = similarities.drop(neuron)

    # Apply threshold if specified
    if threshold is not None:
        similarities = similarities[similarities >= threshold]

    # Sort and take top N
    similarities = similarities.sort_values(ascending=False).head(top_n)

    return pd.DataFrame({
        'neuron_id': similarities.index,
        'similarity': similarities.values
    })


def compare_similarity_dimensions(
    results: dict,
    neuron_pair: tuple[str, str]
) -> pd.DataFrame:
    """
    Compare similarity scores across dimensions for a neuron pair.

    Args:
        results: Output from compute_all_similarities()
        neuron_pair: (neuron_a, neuron_b) tuple

    Returns:
        DataFrame showing similarity in each dimension
    """
    neuron_a, neuron_b = neuron_pair

    comparisons = []
    for dim_name, result in results.items():
        sim_matrix = result['similarity']
        if neuron_a in sim_matrix.index and neuron_b in sim_matrix.columns:
            similarity = sim_matrix.loc[neuron_a, neuron_b]
        else:
            similarity = np.nan

        comparisons.append({
            'dimension': dim_name,
            'similarity': similarity,
            'n_features': result['metadata'].get('n_features', 0)
        })

    return pd.DataFrame(comparisons)


# =============================================================================
# Validation (run with: python -m pub_utils.similarity)
# =============================================================================

if __name__ == "__main__":
    print("Running similarity module validation...")

    # Test 1: Sensory similarity
    print("\n1. Computing sensory similarity...")
    sensory_result = compute_sensory_similarity()
    assert sensory_result['similarity'].shape == (302, 302), "Wrong shape"
    print(f"   Shape: {sensory_result['similarity'].shape}")
    print(f"   Features: {sensory_result['metadata']['n_features']}")

    # Test 2: Process similarity
    print("\n2. Computing process similarity...")
    process_result = compute_process_similarity()
    assert process_result['similarity'].shape == (302, 302), "Wrong shape"
    print(f"   Shape: {process_result['similarity'].shape}")
    print(f"   Features: {process_result['metadata']['n_features']}")

    # Test 3: Verify diagonal = 1.0
    diag = np.diag(sensory_result['similarity'].values)
    assert np.allclose(diag[~np.isnan(diag)], 1.0), "Diagonal should be 1.0"
    print("\n3. Diagonal values verified (all 1.0)")

    # Test 4: Verify symmetry
    sim = sensory_result['similarity'].values
    assert np.allclose(sim, sim.T, equal_nan=True), "Matrix should be symmetric"
    print("4. Symmetry verified")

    # Test 5: Verify range [0, 1]
    non_nan = sim[~np.isnan(sim)]
    assert non_nan.min() >= 0.0 and non_nan.max() <= 1.0, "Values out of range"
    print("5. Value range verified [0, 1]")

    # Test 6: Find similar neurons
    similar = find_similar_neurons(sensory_result['similarity'], 'ASEL', top_n=5)
    print(f"\n6. Top 5 neurons similar to ASEL (sensory):")
    print(similar.to_string(index=False))

    print("\n All validation tests passed!")
