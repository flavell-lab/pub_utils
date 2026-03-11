# Network Metrics: Mathematical Definitions

All metrics are computed by `src/pub_utils/metrics.py`.

Each connectome is represented as a weighted directed graph $G = (V, E, w)$

where $V$ is the set of $N$ neurons,

$E$ is the set of directed edges (excluding self-loops and zero/NaN entries),

and $w: E \to \mathbb{R}^+$ assigns a positive weight to each edge.

---

## 1. Sparsity (`compute_sparsity`)

**Density** is the fraction of possible directed edges that are present:

$$\text{density} = \frac{|E|}{N(N - 1)}$$

The denominator $N(N-1)$ excludes self-loops. **Sparsity** is the complement:

$$\text{sparsity} = 1 - \text{density}$$

Weight distribution statistics are computed over the set of nonzero edge weights $\{w(e) : e \in E\}$:

| Output key      | Definition |
|-----------------|------------|
| `mean_weight`   | $\bar{w} = \frac{1}{\|E\|} \sum_{e \in E} w(e)$ |
| `median_weight` | $\text{median}\bigl(\{w(e) : e \in E\}\bigr)$ |
| `max_weight`    | $\max_{e \in E} w(e)$ |
| `weight_std`    | $\sigma_w = \sqrt{\frac{1}{\|E\|} \sum_{e \in E} (w(e) - \bar{w})^2}$ |

**Code reference:** `metrics.py` lines 78–123. Graph construction in `_to_digraph` (lines 50–71) excludes self-loops (`src == tgt`), zero weights, and NaN entries.

---

## 2. Modularity (`compute_modularity`)

Community detection uses the **Louvain algorithm** (`nx.community.louvain_communities`) on the weighted directed graph, with `resolution=1.0` and `seed=42` for reproducibility.

Given a partition of $V$ into communities $C_1, C_2, \ldots, C_k$, **modularity $Q$** is defined as:

$$Q = \frac{1}{m} \sum_{i,j} \left[ A_{ij} - \frac{s_i^{\text{out}} \, s_j^{\text{in}}}{m} \right] \delta(c_i, c_j)$$

where:

- $A_{ij} = w(i,j)$ if edge $(i,j) \in E$, else $0$
- $m = \sum_{i,j} A_{ij}$ (total edge weight)
- $s_i^{\text{out}} = \sum_j A_{ij}$ (out-strength of node $i$)
- $s_j^{\text{in}} = \sum_i A_{ij}$ (in-strength of node $j$)
- $\delta(c_i, c_j) = 1$ if nodes $i$ and $j$ are in the same community, else $0$

$Q$ ranges from approximately $-0.5$ to $1.0$. Values above ${\sim}0.3$ indicate meaningful community structure.

| Output key            | Definition |
|-----------------------|------------|
| `Q`                   | Modularity score (formula above) |
| `n_communities`       | $k$ = number of communities detected |
| `community_sizes`     | list of $\|C_i\|$ sorted descending |
| `neuron_communities`  | dict mapping each neuron to its community index |

**Code reference:** `metrics.py` lines 126–163. Uses `nx.community.louvain_communities` (line 144) and `nx.community.modularity` (line 148).

---

## 3. Hubness (`compute_hubness`)

Hubness quantifies how unevenly connectivity is distributed across neurons. The primary measure is the **coefficient of variation (CV)** of weighted degree distributions.

**Weighted degrees** for each node $i$:

$$d_i^{\text{in}} = \sum_j w(j, i) \qquad d_i^{\text{out}} = \sum_j w(i, j) \qquad d_i^{\text{total}} = d_i^{\text{in}} + d_i^{\text{out}}$$

**Coefficient of variation:**

$$\text{CV}(\mathbf{x}) = \frac{\sigma(\mathbf{x})}{\mu(\mathbf{x})} = \frac{\sqrt{\frac{1}{N}\sum_{i=1}^{N}(x_i - \bar{x})^2}}{\frac{1}{N}\sum_{i=1}^{N} x_i}$$

A CV near $0$ means uniform degree; a high CV means a few hub neurons dominate the connectivity.

**HITS algorithm** (Kleinberg, 1999) computes authority and hub scores via `nx.hits` with `max_iter=1000` and $L_1$-normalization. For a graph with adjacency matrix $\mathbf{A}$:

$$\mathbf{h} = \mathbf{A}\,\mathbf{a} \qquad \mathbf{a} = \mathbf{A}^\top \mathbf{h}$$

where $\mathbf{h}$ is the hub score vector and $\mathbf{a}$ is the authority score vector. These are the principal eigenvectors of $\mathbf{A}\mathbf{A}^\top$ and $\mathbf{A}^\top\mathbf{A}$, respectively, solved by power iteration.

If convergence fails, a warning is emitted and scores default to empty dicts.

| Output key         | Definition |
|--------------------|------------|
| `in_degree_cv`     | $\text{CV}\bigl(\{d_i^{\text{in}} : i \in V\}\bigr)$ |
| `out_degree_cv`    | $\text{CV}\bigl(\{d_i^{\text{out}} : i \in V\}\bigr)$ |
| `total_degree_cv`  | $\text{CV}\bigl(\{d_i^{\text{total}} : i \in V\}\bigr)$ |
| `top_hubs`         | top-$n$ neurons by $d_i^{\text{total}}$, as (neuron, degree) pairs |
| `hub_scores`       | dict $\{i: h_i\}$ for all $i \in V$ |
| `authority_scores` | dict $\{i: a_i\}$ for all $i \in V$ |

**Code reference:** `metrics.py` lines 166–217. Weighted degrees via `G.in_degree(weight='weight')` (line 183), CV helper at lines 187–191, HITS at line 203.

---

## 4. Recurrence (`compute_recurrence`)

Recurrence captures bidirectional and cyclic connectivity.

### 4a. Reciprocity

**Overall reciprocity** (via `nx.overall_reciprocity`) is the fraction of edges that have a reciprocal partner:

$$r = \frac{\bigl|\{(i,j) \in E : (j,i) \in E\}\bigr|}{|E|}$$

This counts each directed edge independently: if both $(A,B)$ and $(B,A)$ exist, both contribute to the numerator.

### 4b. Reciprocal pairs

A **reciprocal pair** is an unordered pair $\{i, j\}$ where both $(i,j) \in E$ and $(j,i) \in E$. Each pair is scored by the minimum weight of the two directions:

$$\text{pair\_strength}(i, j) = \min\bigl(w(i,j),\; w(j,i)\bigr)$$

The top 20 pairs by pair\_strength are returned.

### 4c. Strongly connected components (SCCs)

A **strongly connected component** is a maximal subset $S \subseteq V$ where every node in $S$ can reach every other node in $S$ via directed paths. Computed via `nx.strongly_connected_components` (Tarjan's algorithm).

**Fraction in non-trivial SCCs** measures what proportion of all neurons participate in cycles:

$$f_{\text{SCC}} = \frac{1}{N} \sum_{\substack{S \in \text{SCCs} \\ |S| > 1}} |S|$$

A non-trivial SCC ($|S| > 1$) guarantees the existence of at least one directed cycle through its members.

| Output key                     | Definition |
|--------------------------------|------------|
| `reciprocity`                  | $r$ (formula above) |
| `n_reciprocal_pairs`           | $\bigl|\bigl\{\{i,j\} : (i,j) \in E \land (j,i) \in E\bigr\}\bigr|$ |
| `top_reciprocal_pairs`         | top 20 pairs by $\min(w(i,j), w(j,i))$ |
| `n_scc`                        | total number of SCCs (including singletons) |
| `scc_sizes`                    | list of $|S|$ for each SCC, sorted descending |
| `fraction_in_nontrivial_scc`   | $f_{\text{SCC}}$ (formula above) |

**Code reference:** `metrics.py` lines 220–268. Reciprocity at line 237, reciprocal pair enumeration at lines 241–247, SCCs at line 253.

---

## 5. Clustering Coefficient (`compute_clustering`)

The **global clustering coefficient** (transitivity) measures the tendency of neurons to form triangles — closed loops of three mutually connected nodes.

A **triplet** is a set of three nodes $(i, j, k)$ connected by at least two edges (when treating the graph as undirected). A triplet is **closed** if all three pairwise edges exist (forming a triangle), and **open** if exactly two of the three edges exist.

$$C = \frac{\text{number of closed triplets}}{\text{number of closed triplets} + \text{number of open triplets}} = \frac{3 \times \Delta}{\tau}$$

where:

- $\Delta$ is the number of triangles in the undirected projection of the graph
- $\tau$ is the number of connected triples (paths of length 2)
- The factor of 3 arises because each triangle contains exactly three closed triplets (one centred on each vertex)

$C$ ranges from $0$ (no triangles; e.g., a bipartite or tree-like graph) to $1$ (every pair of a node's neighbours is also connected).

For directed graphs, edges are treated as undirected for triplet counting: an edge $(i, j)$ and/or $(j, i)$ both contribute a single undirected link between $i$ and $j$.

| Output key                | Definition |
|---------------------------|------------|
| `clustering_coefficient`  | $C$ (formula above) |
| `n_triangles`             | $\Delta$ = number of triangles in the undirected projection |

**Code reference:** `metrics.py` lines 270–293. Uses `nx.transitivity` (line 286) and `nx.triangles` on the undirected projection (line 287).

---

## Graph construction notes

All connectomes are treated as **directed graphs**, including electrical (gap junction) matrices. While gap junctions are biophysically bidirectional, the stored matrices may have slight asymmetries from electron microscopy reconstruction, and uniform DiGraph treatment avoids special-casing. This is why electrical synapse connectomes show reciprocity very close to (or exactly) $1.0$.

Edge inclusion criteria (`_to_digraph`, lines 63–69):
- Self-loops excluded: edges where $i = j$ are skipped
- Zero weights excluded: $w(i,j) = 0$ treated as no connection
- NaN excluded: missing data treated as no connection
