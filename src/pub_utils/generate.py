"""
Graph Generators for Metric Showcasing

For each of the 5 network metrics (sparsity/density, modularity, hubness,
reciprocity, clustering), provides a function that generates a small directed
graph with a distinctive value of that metric, then visualizes it as a
2-panel figure (adjacency heatmap + network diagram).

Uses graspologic for simulation/heatmap and networkx for layout/drawing.
"""

import warnings
from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from graspologic.plot import heatmap as gs_heatmap
from graspologic.simulations import er_nm, sbm


# =============================================================================
# Visualization helper
# =============================================================================

def _visualize(adj, title, seed=42, save_path=None):
    """
    Render a 2-panel figure: graspologic heatmap (left) + networkx diagram (right).

    Args:
        adj: N×N numpy adjacency matrix
        title: Figure title
        seed: Random seed for spring layout
        save_path: Optional path to save figure

    Returns:
        matplotlib Figure
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    # Left: graspologic heatmap with gridlines
    N = adj.shape[0]
    gs_heatmap(adj, ax=ax1, title="Adjacency matrix", cbar=False)
    ax1.set_xticks(np.arange(-0.5, N, 1), minor=True)
    ax1.set_yticks(np.arange(-0.5, N, 1), minor=True)
    ax1.grid(which="minor", color="white", linewidth=0.8)
    ax1.tick_params(which="minor", size=0)

    # Right: networkx spring layout
    G = nx.from_numpy_array(adj, create_using=nx.DiGraph)
    pos = nx.spring_layout(G, seed=seed)
    nx.draw(
        G, pos, ax=ax2, with_labels=True,
        node_color="steelblue", node_size=500,
        font_color="white", font_weight="bold",
        edge_color="gray", arrowsize=12,
        connectionstyle="arc3,rad=0.1",
    )
    ax2.set_title("Network diagram")

    fig.suptitle(title, fontsize=13, fontweight="bold")
    fig.tight_layout()

    if save_path:
        out = Path(save_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=150, bbox_inches="tight")

    return fig


# =============================================================================
# Generator functions
# =============================================================================

def generate_sparse(N=5, E=10, seed=42, save_path=None) -> dict:
    """
    Generate a sparse Erdos-Renyi directed graph.

    Showcases density = E / (N*(N-1)).

    Args:
        N: Number of nodes
        E: Number of edges
        seed: Random seed
        save_path: Optional path to save the figure

    Returns:
        dict with adjacency, metric_name, metric_value, figure
    """
    np.random.seed(seed)
    adj = er_nm(N, E, directed=True)

    n_edges = int(adj.sum())
    max_possible = N * (N - 1)
    density = n_edges / max_possible if max_possible > 0 else 0.0

    title = f"Sparse graph: density = {density:.3f}  (N={N}, E={n_edges})"
    fig = _visualize(adj, title, seed=seed, save_path=save_path)

    return {
        "adjacency": adj,
        "metric_name": "density",
        "metric_value": density,
        "figure": fig,
    }


def generate_modular(N=6, E=12, seed=42, save_path=None) -> dict:
    """
    Generate a stochastic block model graph with clear community structure.

    Partitions N nodes into 2 blocks with high intra-block and low inter-block
    connection probabilities, tuned so expected edge count ~= E.

    Showcases modularity Q.

    Args:
        N: Number of nodes (split into 2 blocks)
        E: Number of edges (approximate)
        seed: Random seed
        save_path: Optional path to save the figure

    Returns:
        dict with adjacency, metric_name, metric_value, figure
    """
    np.random.seed(seed)

    # Split into 2 blocks
    n1 = N // 2
    n2 = N - n1
    n_per_block = [n1, n2]

    # For small N, manually build block-diagonal adjacency to guarantee high Q
    intra_pairs = n1 * (n1 - 1) + n2 * (n2 - 1)
    if intra_pairs >= E:
        # All edges can fit intra-block: use p_inter=0
        p_intra = min(E / intra_pairs, 1.0)
        P = np.array([[p_intra, 0.0], [0.0, p_intra]])
        adj = sbm(n_per_block, P, directed=True)
    else:
        # Not enough intra-block capacity: fill intra-block completely,
        # add remaining edges inter-block
        adj = np.zeros((N, N), dtype=int)
        block1 = list(range(n1))
        block2 = list(range(n1, N))
        # Fill all intra-block edges
        for blk in [block1, block2]:
            for i in blk:
                for j in blk:
                    if i != j:
                        adj[i, j] = 1
        # Add inter-block edges up to E
        n_placed = int(adj.sum())
        inter_candidates = [(i, j) for i in block1 for j in block2] + \
                           [(i, j) for i in block2 for j in block1]
        np.random.shuffle(inter_candidates)
        for i, j in inter_candidates:
            if n_placed >= E:
                break
            adj[i, j] = 1
            n_placed += 1

    # Compute modularity via networkx
    G = nx.from_numpy_array(adj, create_using=nx.DiGraph)
    communities = nx.community.louvain_communities(G, seed=seed)
    Q = nx.community.modularity(G, communities)

    title = f"Modular graph: Q = {Q:.3f}  (N={N}, {len(communities)} communities)"
    fig = _visualize(adj, title, seed=seed, save_path=save_path)

    return {
        "adjacency": adj,
        "metric_name": "modularity",
        "metric_value": Q,
        "figure": fig,
    }


def generate_hub(N=5, E=10, seed=42, save_path=None) -> dict:
    """
    Generate a hub-dominated directed graph (star-like topology).

    One hub node connects to/from all others; remaining edges are distributed
    randomly among non-hub nodes.

    Showcases total degree CV (coefficient of variation).

    Args:
        N: Number of nodes
        E: Number of edges
        seed: Random seed
        save_path: Optional path to save the figure

    Returns:
        dict with adjacency, metric_name, metric_value, figure
    """
    np.random.seed(seed)
    adj = np.zeros((N, N), dtype=int)

    hub = 0
    edges_placed = 0

    # Hub -> all others (outgoing)
    for j in range(1, N):
        if edges_placed >= E:
            break
        adj[hub, j] = 1
        edges_placed += 1

    # All others -> hub (incoming)
    for i in range(1, N):
        if edges_placed >= E:
            break
        adj[i, hub] = 1
        edges_placed += 1

    # Fill remaining edges randomly among non-hub nodes
    if edges_placed < E:
        non_hub = list(range(1, N))
        candidates = [(i, j) for i in non_hub for j in non_hub if i != j and adj[i, j] == 0]
        np.random.shuffle(candidates)
        for i, j in candidates:
            if edges_placed >= E:
                break
            adj[i, j] = 1
            edges_placed += 1

    # Compute degree CV
    G = nx.from_numpy_array(adj, create_using=nx.DiGraph)
    in_deg = np.array([d for _, d in G.in_degree()])
    out_deg = np.array([d for _, d in G.out_degree()])
    total_deg = in_deg + out_deg
    mean_deg = np.mean(total_deg)
    cv = float(np.std(total_deg) / mean_deg) if mean_deg > 0 else 0.0

    title = f"Hub graph: degree CV = {cv:.3f}  (N={N}, E={int(adj.sum())})"
    fig = _visualize(adj, title, seed=seed, save_path=save_path)

    return {
        "adjacency": adj,
        "metric_name": "total_degree_cv",
        "metric_value": cv,
        "figure": fig,
    }


def generate_reciprocal(N=5, E=10, seed=42, save_path=None) -> dict:
    """
    Generate a directed graph with high reciprocity.

    Creates E//2 undirected edges, then symmetrizes them to form reciprocal
    directed pairs. Adds one extra edge if E is odd.

    Showcases reciprocity ~= 1.0.

    Args:
        N: Number of nodes
        E: Number of edges (directed count after symmetrization)
        seed: Random seed
        save_path: Optional path to save the figure

    Returns:
        dict with adjacency, metric_name, metric_value, figure
    """
    np.random.seed(seed)

    # Generate undirected edges (E//2 of them), then symmetrize
    n_undirected = E // 2
    max_undirected = N * (N - 1) // 2
    n_undirected = min(n_undirected, max_undirected)

    adj_undirected = er_nm(N, n_undirected, directed=False)

    # Symmetrize to create directed reciprocal pairs
    adj = np.maximum(adj_undirected, adj_undirected.T).astype(int)

    # Add one extra random edge if E is odd
    if E % 2 == 1:
        zeros = list(zip(*np.where((adj == 0) & (np.eye(N) == 0))))
        if zeros:
            idx = np.random.randint(len(zeros))
            i, j = zeros[idx]
            adj[i, j] = 1

    # Compute reciprocity
    G = nx.from_numpy_array(adj, create_using=nx.DiGraph)
    n_edges = G.number_of_edges()
    reciprocity = nx.overall_reciprocity(G) if n_edges > 0 else 0.0

    title = f"Reciprocal graph: reciprocity = {reciprocity:.3f}  (N={N}, E={n_edges})"
    fig = _visualize(adj, title, seed=seed, save_path=save_path)

    return {
        "adjacency": adj,
        "metric_name": "reciprocity",
        "metric_value": reciprocity,
        "figure": fig,
    }


def generate_clustered(N=5, E=10, seed=42, save_path=None) -> dict:
    """
    Generate a directed graph with high clustering coefficient.

    Builds triangles greedily: picks random triples of nodes and adds all 6
    directed edges (both directions for each of 3 undirected edges in the
    triangle). For small N, starts from a clique on a subset of nodes.

    Showcases global clustering coefficient (transitivity).

    Args:
        N: Number of nodes
        E: Number of edges (approximate)
        seed: Random seed
        save_path: Optional path to save the figure

    Returns:
        dict with adjacency, metric_name, metric_value, figure
    """
    np.random.seed(seed)
    adj = np.zeros((N, N), dtype=int)
    np.fill_diagonal(adj, 0)

    nodes = list(range(N))
    edges_placed = int(adj.sum())

    # Greedily add triangles
    max_attempts = 500
    attempt = 0
    while edges_placed < E and attempt < max_attempts:
        attempt += 1
        triple = np.random.choice(nodes, size=3, replace=False)
        a, b, c = triple

        # Add all 6 directed edges for this triangle
        for i, j in [(a, b), (b, a), (b, c), (c, b), (a, c), (c, a)]:
            if adj[i, j] == 0 and edges_placed < E:
                adj[i, j] = 1
                edges_placed += 1

    # Compute clustering coefficient (transitivity)
    G = nx.from_numpy_array(adj, create_using=nx.DiGraph)
    cc = nx.transitivity(G)

    title = f"Clustered graph: CC = {cc:.3f}  (N={N}, E={int(adj.sum())})"
    fig = _visualize(adj, title, seed=seed, save_path=save_path)

    return {
        "adjacency": adj,
        "metric_name": "clustering_coefficient",
        "metric_value": cc,
        "figure": fig,
    }
