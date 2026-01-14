import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import seaborn as sns
from .constants import AllHermNeurons, AllHermNeuronBlocks

def plot_connectome_matrix(plot_df, title="", colormap_name='hot', colorbar_label='# Unique Ligand-Receptor Pairs', colormap_thresh=None, show_blocks=True):
    """
    Plot a connectome matrix as a heatmap.

    Args:
        plot_df: DataFrame with connectome data
        title: Plot title
        colormap_name: Name of matplotlib colormap to use
        colorbar_label: Label for the colorbar
        colormap_thresh: If max value exceeds this threshold, use 97th percentile as colormap max
        show_blocks: If True and neuron order matches AllHermNeurons, show block dividers and labels
    """
    # Determine max value and create adaptive colormap from continuous colormap
    actual_max = int(plot_df.max().max())

    # Check if we need to clip the colormap
    if colormap_thresh is not None and actual_max > colormap_thresh:
        # Compute 90th percentile of non-zero values
        values = plot_df.values.flatten()
        values = values[~np.isnan(values)]
        values = values[values > 0]
        max_val = int(np.percentile(values, 97))
        print(f"Actual max value: {actual_max}, using 97th percentile ({max_val}) as colormap max")
        colormap_clipped = True
    else:
        max_val = actual_max
        colormap_clipped = False

    num_colors = max_val + 1  # Include 0
    
    # Get the base colormap
    base_cmap = plt.get_cmap(colormap_name)
    
    # Sample colors evenly from the colormap
    if num_colors == 1:
        color_indices = [0]
    else:
        # Sample from 0 to 0.9 of the colormap range for better contrast
        color_indices = np.linspace(0, 0.9, num_colors)
    
    colors = [base_cmap(idx) for idx in color_indices]
    
    # Create discrete colormap
    cmap = mcolors.ListedColormap(colors)
    
    # --- Modification: Set NaN values to grey ---
    cmap.set_bad("grey") 
    # --------------------------------------------

    bounds = list(range(max_val + 2))  # [0, 1, 2, ..., max_val+1]
    norm = mcolors.BoundaryNorm(bounds, cmap.N, clip=True)

    # Clip data to max_val if colormap is clipped (so values > max_val render as max color)
    if colormap_clipped:
        plot_df = plot_df.clip(upper=max_val)
    
    # Create tick positions and labels for colorbar
    tick_positions = [i + 0.5 for i in range(max_val + 1)]
    tick_labels = [str(i) for i in range(max_val + 1)]

    # If colormap was clipped, indicate that top value is "X to Y"
    if colormap_clipped:
        tick_labels[-1] = f"{max_val} to {actual_max}"

    # 4. Plotting
    fig, ax = plt.subplots(figsize=(20, 22))
    
    # Set aspect ratio to 1
    ax.set_aspect('equal')

    sns.heatmap(
        plot_df, 
        fmt="", 
        cmap=cmap,
        norm=norm,
        cbar=True,
        cbar_kws={
            'orientation': 'horizontal',
            'shrink': 0.4,
            'pad': 0.05,
            'ticks': tick_positions
        },
        xticklabels=True, 
        yticklabels=True,
        linewidths=0.05,  # Thinner gridlines
        linecolor=(0.5, 0.5, 0.5, 0.2),  # Grey with transparency (RGBA)
        ax=ax
    )
    
    # --- Customizations ---
    # Colorbar modifications (Boxed and Bold Label)
    cbar = ax.collections[0].colorbar
    cbar.outline.set_visible(True)
    cbar.outline.set_linewidth(1)
    cbar.outline.set_edgecolor('black')
    
    cbar.set_ticklabels(tick_labels)
    
    # Make colorbar label larger and bolded
    cbar.set_label(colorbar_label, size=20, weight='bold', labelpad=5)
    
    # Titles and Labels
    plt.title(title, fontsize=30, pad=36, fontweight='bold')
    plt.xticks(rotation=90, fontsize=5)
    plt.yticks(rotation=0, fontsize=5)
    plt.xlabel(f'{len(plot_df.columns)} Recipient Neurons', fontsize=20, fontweight='bold')
    plt.ylabel(f'{len(plot_df.index)} Source Neurons', fontsize=20, fontweight='bold')

    # Check if neuron order matches AllHermNeurons and add block dividers/labels
    neuron_order = list(plot_df.index)
    if show_blocks and neuron_order == AllHermNeurons:
        n_neurons = len(neuron_order)

        # Draw white lines at block boundaries
        for block_name, start_idx, end_idx in AllHermNeuronBlocks[:-1]:  # Skip last block (no line after it)
            boundary = end_idx + 1  # Line after the last neuron of this block
            # Horizontal line (for rows/source neurons)
            ax.axhline(y=boundary, color='white', linewidth=2, zorder=10)
            # Vertical line (for columns/recipient neurons)
            ax.axvline(x=boundary, color='white', linewidth=2, zorder=10)

        # Add block labels on the sides
        for block_name, start_idx, end_idx in AllHermNeuronBlocks:
            mid_idx = (start_idx + end_idx) / 2 + 0.5  # Center of the block

            # Label on the right side (for rows)
            ax.text(n_neurons + 1, mid_idx, block_name,
                    va='center', ha='left', fontsize=20, fontweight='bold',
                    rotation=270)

            # Label on the top (for columns)
            ax.text(mid_idx, -1, block_name,
                    va='bottom', ha='center', fontsize=20, fontweight='bold',
                    rotation=0)

    plt.tight_layout()
    plt.show()
    
    # Return the figure
    return fig


import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

def plot_reciprocal_network(reciprocal_pairs, figsize=(10, 10), title='', colormap_name='hot'):
    """
    Create a network diagram showing reciprocal connections.
    Arrow thickness and color indicate connection strength (v12, v21 values).
    Uses the same colormap as plot_connectome_matrix.
    """
    G = nx.DiGraph()
    
    # Add edges with weights
    for n1, n2, v12, v21 in reciprocal_pairs:
        G.add_edge(n1, n2, weight=v12)
        G.add_edge(n2, n1, weight=v21)
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Use Kamada-Kawai layout for better spacing
    pos = nx.kamada_kawai_layout(G)
    
    # Draw nodes - black and white style
    node_sizes = [400 + 120 * G.degree(n) for n in G.nodes()]
    nx.draw_networkx_nodes(G, pos, node_size=node_sizes, 
                            node_color='white', 
                            edgecolors='black', 
                            linewidths=2, ax=ax)
    
    # Draw labels - larger text
    nx.draw_networkx_labels(G, pos, font_size=9, font_weight='bold', 
                            font_color='black', ax=ax)
    
    # Get weights for edge styling
    weights = np.array([G[u][v]['weight'] for u, v in G.edges()])
    max_weight = int(weights.max())
    min_weight = int(weights.min())
    
    # Create discrete colormap matching plot_connectome_matrix
    base_cmap = plt.get_cmap(colormap_name)
    num_colors = max_weight + 1
    
    if num_colors == 1:
        color_indices = [0.9]  # Use a visible color for single value
    else:
        # Sample from 0 to 0.9 of the colormap range for better contrast
        color_indices = np.linspace(0, 0.9, num_colors)
    
    colors = [base_cmap(idx) for idx in color_indices]
    cmap = mcolors.ListedColormap(colors)
    bounds = list(range(max_weight + 2))
    norm = mcolors.BoundaryNorm(bounds, cmap.N)
    
    # Draw edges with arrows - colored by weight using hot colormap
    for (u, v), weight in zip(G.edges(), weights):
        edge_color = cmap(norm(weight))
        nx.draw_networkx_edges(G, pos,
                                edgelist=[(u, v)],
                                width=1.5 + 2 * (weight / max_weight),
                                alpha=0.9,
                                edge_color=[edge_color],
                                arrows=True,
                                arrowsize=20,
                                arrowstyle='-|>',
                                connectionstyle="arc3,rad=0.15",
                                ax=ax,
                                node_size=node_sizes[0])  # For proper arrow positioning
    
    # Create discrete colorbar matching plot_connectome_matrix style
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    
    # Position colorbar horizontally at bottom
    cbar = plt.colorbar(sm, ax=ax, orientation='horizontal', 
                        shrink=0.3, pad=0.05)
    
    # Style colorbar to match plot_connectome_matrix
    cbar.outline.set_visible(True)
    cbar.outline.set_linewidth(1)
    cbar.outline.set_edgecolor('black')
    
    # Set discrete ticks at center of each color band
    tick_positions = [i + 0.5 for i in range(max_weight + 1)]
    tick_labels = [str(i) for i in range(max_weight + 1)]
    cbar.set_ticks(tick_positions)
    cbar.set_ticklabels(tick_labels)
    cbar.set_label('# Unique Ligand-Receptor Pairs', size=14, weight='bold', labelpad=10)
    
    ax.set_title(title, 
                 fontsize=20, fontweight='bold', pad=20)
    ax.axis('off')
    
    plt.tight_layout()
    return fig, ax, G