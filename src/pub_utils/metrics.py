"""
Network Metrics Module

Computes graph-theoretic network metrics for C. elegans connectomes:
1. Sparsity - Edge density and weight distribution statistics
2. Modularity - Community structure via Louvain algorithm
3. Hubness - Degree distribution asymmetry and hub identification
4. Recurrence - Reciprocal connections and strongly connected components
5. Clustering - Global clustering coefficient (transitivity)

All functions accept either a pd.DataFrame (302×302 adjacency matrix) or a
NeuronInteraction object. Connectomes are treated uniformly as directed graphs.
"""

import json
import warnings
from pathlib import Path
from typing import Optional, Union

import networkx as nx
import numpy as np
import pandas as pd

from .core import NeuronInteraction

# Base path points to repository root
_BASE_PATH = Path(__file__).parent.parent.parent

# Module-level cache
_ASSETS_CACHE = None


# =============================================================================
# Internal Helpers
# =============================================================================

def _load_assets() -> dict:
    """Load and cache assets.json."""
    global _ASSETS_CACHE
    if _ASSETS_CACHE is None:
        with open(_BASE_PATH / "assets.json") as f:
            _ASSETS_CACHE = json.load(f)
    return _ASSETS_CACHE


def _get_path(asset_path: str) -> Path:
    """Convert relative asset path to absolute path."""
    return _BASE_PATH / asset_path


def _to_digraph(matrix: Union[pd.DataFrame, NeuronInteraction]) -> nx.DiGraph:
    """
    Convert a DataFrame or NeuronInteraction to a weighted nx.DiGraph.

    Self-loops are excluded. Zero-weight edges are excluded.
    NaN values are treated as absent edges.
    """
    if isinstance(matrix, NeuronInteraction):
        matrix = matrix.data

    G = nx.DiGraph()
    G.add_nodes_from(matrix.index)

    for src in matrix.index:
        for tgt in matrix.columns:
            if src == tgt:
                continue
            w = matrix.loc[src, tgt]
            if pd.notna(w) and w != 0:
                G.add_edge(src, tgt, weight=float(w))

    return G


# =============================================================================
# Metric Functions
# =============================================================================

def compute_sparsity(matrix: Union[pd.DataFrame, NeuronInteraction]) -> dict:
    """
    Compute sparsity and weight distribution statistics.

    Excludes self-loops: density = |E| / (N × (N-1)).

    Args:
        matrix: 302×302 weighted adjacency matrix or NeuronInteraction

    Returns:
        dict with keys: sparsity, density, n_nodes, n_edges,
        mean_weight, median_weight, max_weight, weight_std
    """
    G = _to_digraph(matrix)
    n = G.number_of_nodes()
    n_edges = G.number_of_edges()
    max_possible = n * (n - 1)

    density = n_edges / max_possible if max_possible > 0 else 0.0
    sparsity = 1.0 - density

    weights = [d['weight'] for _, _, d in G.edges(data=True)]

    if weights:
        w_arr = np.array(weights)
        stats = {
            'mean_weight': float(np.mean(w_arr)),
            'median_weight': float(np.median(w_arr)),
            'max_weight': float(np.max(w_arr)),
            'weight_std': float(np.std(w_arr)),
        }
    else:
        stats = {
            'mean_weight': 0.0,
            'median_weight': 0.0,
            'max_weight': 0.0,
            'weight_std': 0.0,
        }

    return {
        'sparsity': sparsity,
        'density': density,
        'n_nodes': n,
        'n_edges': n_edges,
        **stats,
    }


def compute_modularity(
    matrix: Union[pd.DataFrame, NeuronInteraction],
    resolution: float = 1.0,
    seed: int = 42,
) -> dict:
    """
    Detect community structure via Louvain algorithm and compute modularity Q.

    Args:
        matrix: 302×302 weighted adjacency matrix or NeuronInteraction
        resolution: Louvain resolution parameter (higher = more communities)
        seed: Random seed for reproducibility

    Returns:
        dict with keys: Q, n_communities, community_sizes, neuron_communities
    """
    G = _to_digraph(matrix)

    communities = nx.community.louvain_communities(
        G, weight='weight', resolution=resolution, seed=seed
    )

    Q = nx.community.modularity(G, communities, weight='weight')

    # Build neuron → community mapping
    neuron_communities = {}
    community_sizes = []
    for i, comm in enumerate(communities):
        community_sizes.append(len(comm))
        for neuron in comm:
            neuron_communities[neuron] = i

    return {
        'Q': Q,
        'n_communities': len(communities),
        'community_sizes': sorted(community_sizes, reverse=True),
        'neuron_communities': neuron_communities,
    }


def compute_hubness(
    matrix: Union[pd.DataFrame, NeuronInteraction],
    top_n: int = 10,
) -> dict:
    """
    Quantify hub-dominated topology via degree distribution CV and HITS scores.

    Args:
        matrix: 302×302 weighted adjacency matrix or NeuronInteraction
        top_n: Number of top hub neurons to return

    Returns:
        dict with keys: in_degree_cv, out_degree_cv, total_degree_cv,
        top_hubs (list of (neuron, total_degree)), hub_scores, authority_scores
    """
    G = _to_digraph(matrix)

    in_deg = np.array([d for _, d in G.in_degree(weight='weight')])
    out_deg = np.array([d for _, d in G.out_degree(weight='weight')])
    total_deg = in_deg + out_deg

    def _cv(arr):
        mean = np.mean(arr)
        if mean == 0:
            return 0.0
        return float(np.std(arr) / mean)

    # Top hub neurons by total weighted degree
    nodes = list(G.nodes())
    node_total = [(nodes[i], float(total_deg[i])) for i in range(len(nodes))]
    node_total.sort(key=lambda x: x[1], reverse=True)
    top_hubs = node_total[:top_n]

    # HITS algorithm
    hub_scores = {}
    authority_scores = {}
    try:
        hub_scores, authority_scores = nx.hits(G, max_iter=1000, normalized=True)
    except nx.PowerIterationFailedConvergence:
        warnings.warn(
            "HITS algorithm did not converge; hub_scores and authority_scores "
            "will be empty dicts"
        )

    return {
        'in_degree_cv': _cv(in_deg),
        'out_degree_cv': _cv(out_deg),
        'total_degree_cv': _cv(total_deg),
        'top_hubs': top_hubs,
        'hub_scores': hub_scores,
        'authority_scores': authority_scores,
    }


def compute_recurrence(matrix: Union[pd.DataFrame, NeuronInteraction]) -> dict:
    """
    Quantify reciprocal connectivity and cycle participation.

    Args:
        matrix: 302×302 weighted adjacency matrix or NeuronInteraction

    Returns:
        dict with keys: reciprocity, n_reciprocal_pairs, top_reciprocal_pairs,
        n_scc, scc_sizes, fraction_in_nontrivial_scc
    """
    G = _to_digraph(matrix)

    # Reciprocity
    if G.number_of_edges() == 0:
        reciprocity = 0.0
    else:
        reciprocity = nx.overall_reciprocity(G)

    # Enumerate reciprocal pairs with min(w_AB, w_BA)
    reciprocal_pairs = []
    seen = set()
    for u, v, d in G.edges(data=True):
        if G.has_edge(v, u) and (v, u) not in seen:
            seen.add((u, v))
            w_uv = d['weight']
            w_vu = G[v][u]['weight']
            reciprocal_pairs.append((u, v, min(w_uv, w_vu)))

    reciprocal_pairs.sort(key=lambda x: x[2], reverse=True)
    top_reciprocal = reciprocal_pairs[:20]

    # Strongly connected components
    sccs = list(nx.strongly_connected_components(G))
    scc_sizes = sorted([len(s) for s in sccs], reverse=True)

    # Fraction of nodes in non-trivial SCCs (size > 1)
    n_in_nontrivial = sum(len(s) for s in sccs if len(s) > 1)
    n_total = G.number_of_nodes()
    frac_nontrivial = n_in_nontrivial / n_total if n_total > 0 else 0.0

    return {
        'reciprocity': reciprocity,
        'n_reciprocal_pairs': len(reciprocal_pairs),
        'top_reciprocal_pairs': top_reciprocal,
        'n_scc': len(sccs),
        'scc_sizes': scc_sizes,
        'fraction_in_nontrivial_scc': frac_nontrivial,
    }


def compute_clustering(matrix: Union[pd.DataFrame, NeuronInteraction]) -> dict:
    """
    Compute the global clustering coefficient (transitivity).

    Transitivity = (number of closed triplets) / (number of closed triplets +
    number of open triplets), equivalent to 3 * triangles / connected triples.

    For directed graphs, edges are treated as undirected for triplet counting.

    Args:
        matrix: 302×302 weighted adjacency matrix or NeuronInteraction

    Returns:
        dict with keys: clustering_coefficient, n_triangles
    """
    G = _to_digraph(matrix)

    clustering_coefficient = nx.transitivity(G)
    n_triangles = sum(nx.triangles(G.to_undirected()).values()) // 3

    return {
        'clustering_coefficient': clustering_coefficient,
        'n_triangles': n_triangles,
    }


# =============================================================================
# Convenience Functions
# =============================================================================

def compute_metrics(matrix: Union[pd.DataFrame, NeuronInteraction]) -> dict:
    """
    Compute all four metric categories for a single connectome.

    Returns:
        dict with keys: 'sparsity', 'modularity', 'hubness', 'recurrence'
    """
    return {
        'sparsity': compute_sparsity(matrix),
        'modularity': compute_modularity(matrix),
        'hubness': compute_hubness(matrix),
        'recurrence': compute_recurrence(matrix),
        'clustering': compute_clustering(matrix),
    }


def compute_all_metrics(save_path: Optional[str] = None) -> pd.DataFrame:
    """
    Compute scalar network metrics for all 45 preassembled connectomes.

    Iterates over structural (chemical/electrical) and molecular (neuropeptide/
    monoamine) connectomes listed in assets.json.

    Args:
        save_path: Optional path to save the summary CSV

    Returns:
        DataFrame with one row per connectome and scalar metric columns:
        connectome_type, dataset, n_edges, density, sparsity, mean_weight,
        max_weight, Q, n_communities, in_degree_cv, out_degree_cv,
        total_degree_cv, reciprocity, n_reciprocal_pairs,
        fraction_in_nontrivial_scc, clustering_coefficient, n_triangles
    """
    assets = _load_assets()
    rows = []

    # --- Structural connectomes ---
    for synapse_type in ['chemical_synapse', 'electrical_synapse']:
        datasets = assets['structural_connectomes']['preassembled'].get(synapse_type, {})
        for dataset_name, csv_path in datasets.items():
            label = f"{synapse_type}_{dataset_name}"
            print(f"  Processing {label}...")
            df = pd.read_csv(_get_path(csv_path), index_col=0)
            rows.append(_metrics_row(label, synapse_type, dataset_name, df))

    # --- Molecular connectomes ---
    for mol_type in ['neuropeptide', 'monoamine']:
        datasets = assets['molecular_connectomes']['preassembled'].get(mol_type, {})
        for dataset_name, csv_path in datasets.items():
            label = f"{mol_type}_{dataset_name}"
            print(f"  Processing {label}...")
            df = pd.read_csv(_get_path(csv_path), index_col=0)
            rows.append(_metrics_row(label, mol_type, dataset_name, df))

    summary = pd.DataFrame(rows)

    if save_path:
        out = Path(save_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        summary.to_csv(out, index=False)
        print(f"Saved summary to {out}")

    return summary


def _metrics_row(label: str, conn_type: str, dataset: str, df: pd.DataFrame) -> dict:
    """Build a single summary row for compute_all_metrics."""
    sp = compute_sparsity(df)
    mod = compute_modularity(df)
    hub = compute_hubness(df)
    rec = compute_recurrence(df)
    cl = compute_clustering(df)

    return {
        'label': label,
        'connectome_type': conn_type,
        'dataset': dataset,
        'n_nodes': sp['n_nodes'],
        'n_edges': sp['n_edges'],
        'density': sp['density'],
        'sparsity': sp['sparsity'],
        'mean_weight': sp['mean_weight'],
        'max_weight': sp['max_weight'],
        'Q': mod['Q'],
        'n_communities': mod['n_communities'],
        'in_degree_cv': hub['in_degree_cv'],
        'out_degree_cv': hub['out_degree_cv'],
        'total_degree_cv': hub['total_degree_cv'],
        'reciprocity': rec['reciprocity'],
        'n_reciprocal_pairs': rec['n_reciprocal_pairs'],
        'fraction_in_nontrivial_scc': rec['fraction_in_nontrivial_scc'],
        'clustering_coefficient': cl['clustering_coefficient'],
        'n_triangles': cl['n_triangles'],
    }


# =============================================================================
# Validation (run with: python -m pub_utils.metrics)
# =============================================================================

if __name__ == "__main__":
    print("Running metrics module validation...")

    # Load a well-characterized structural connectome
    assets = _load_assets()
    cook_path = _get_path(
        assets['structural_connectomes']['preassembled']['chemical_synapse']['Cook2019']
    )
    cook = pd.read_csv(cook_path, index_col=0)
    print(f"\nLoaded Cook2019 chemical connectome: {cook.shape}")

    # Test 1: Sparsity
    print("\n1. Sparsity metrics:")
    sp = compute_sparsity(cook)
    print(f"   Density:      {sp['density']:.4f}")
    print(f"   Sparsity:     {sp['sparsity']:.4f}")
    print(f"   Edges:        {sp['n_edges']}")
    print(f"   Mean weight:  {sp['mean_weight']:.2f}")
    print(f"   Max weight:   {sp['max_weight']:.0f}")
    assert 0 < sp['density'] < 1, "Density out of range"
    assert sp['sparsity'] == 1.0 - sp['density'], "Sparsity/density mismatch"

    # Test 2: Modularity
    print("\n2. Modularity metrics:")
    mod = compute_modularity(cook)
    print(f"   Q:              {mod['Q']:.4f}")
    print(f"   Communities:    {mod['n_communities']}")
    print(f"   Largest comm:   {mod['community_sizes'][0]}")
    assert -0.5 <= mod['Q'] <= 1.0, "Q out of theoretical range"
    assert sum(mod['community_sizes']) == cook.shape[0], "Community sizes don't sum to N"

    # Test 3: Hubness
    print("\n3. Hubness metrics:")
    hub = compute_hubness(cook)
    print(f"   In-degree CV:   {hub['in_degree_cv']:.3f}")
    print(f"   Out-degree CV:  {hub['out_degree_cv']:.3f}")
    print(f"   Total-deg CV:   {hub['total_degree_cv']:.3f}")
    print(f"   Top 5 hubs:     {[(n, f'{d:.0f}') for n, d in hub['top_hubs'][:5]]}")
    assert hub['total_degree_cv'] > 0, "CV should be positive for a real network"

    # Test 4: Recurrence
    print("\n4. Recurrence metrics:")
    rec = compute_recurrence(cook)
    print(f"   Reciprocity:    {rec['reciprocity']:.4f}")
    print(f"   Recip. pairs:   {rec['n_reciprocal_pairs']}")
    print(f"   SCCs:           {rec['n_scc']}")
    print(f"   Largest SCC:    {rec['scc_sizes'][0]}")
    print(f"   Frac nontrivial SCC: {rec['fraction_in_nontrivial_scc']:.3f}")
    assert 0 <= rec['reciprocity'] <= 1, "Reciprocity out of range"

    # Test 5: Combined compute_metrics
    print("\n5. compute_metrics (combined):")
    all_m = compute_metrics(cook)
    assert set(all_m.keys()) == {'sparsity', 'modularity', 'hubness', 'recurrence'}
    print("   All four metric dicts returned correctly")

    # Test 6: NeuronInteraction compatibility
    print("\n6. NeuronInteraction compatibility:")
    ni = NeuronInteraction(cook)
    ni_sp = compute_sparsity(ni)
    assert ni_sp['n_edges'] == sp['n_edges'], "NeuronInteraction should give same result"
    print("   NeuronInteraction input matches DataFrame input")

    print("\n All validation tests passed!")
