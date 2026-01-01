import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors

def plot_connectome_matrix(df, title="", colormap_name='hot'):
    plot_df = df.copy()
    
    # 1. Handle the 'Row' column and clean neuron names
    if 'Row' in plot_df.columns:
        plot_df['Row'] = plot_df['Row'].astype(str).str.strip()
        plot_df = plot_df.set_index('Row')
    
    # 2. Synchronize Ticks: Sorted alphabetically
    plot_df.index = plot_df.index.astype(str)
    plot_df.columns = plot_df.columns.astype(str)
    all_neurons = sorted(list(set(plot_df.index) | set(plot_df.columns)))
    plot_df = plot_df.reindex(index=all_neurons, columns=all_neurons).fillna(0)

    # 3. Determine max value and create adaptive colormap from continuous colormap
    max_val = int(plot_df.max().max())
    num_colors = max_val + 1  # Include 0
    
    # Get the base colormap
    base_cmap = plt.get_cmap(colormap_name)
    
    # Sample colors evenly from the colormap
    # Use linspace to get evenly spaced values, avoiding the extreme end (1.0) 
    # which can be too light in some colormaps
    if num_colors == 1:
        color_indices = [0]
    else:
        # Sample from 0 to 0.9 of the colormap range for better contrast
        color_indices = np.linspace(0, 0.9, num_colors)
    
    colors = [base_cmap(idx) for idx in color_indices]
    
    # Create discrete colormap
    cmap = mcolors.ListedColormap(colors)
    bounds = list(range(max_val + 2))  # [0, 1, 2, ..., max_val+1]
    norm = mcolors.BoundaryNorm(bounds, cmap.N)
    
    # Create tick positions and labels for colorbar
    tick_positions = [i + 0.5 for i in range(max_val + 1)]
    tick_labels = [str(i) for i in range(max_val + 1)]

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
            'shrink': 0.3,
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
    cbar.set_label('# Ligand-Receptor Pairs', size=16, weight='bold', labelpad=5)
    
    # Titles and Labels
    plt.title(title, fontsize=20, pad=10, fontweight='bold')
    plt.xticks(rotation=90, fontsize=5)
    plt.yticks(rotation=0, fontsize=5)
    plt.xlabel(f'{len(plot_df.columns)} Recipient Neurons', fontsize=16, fontweight='bold')
    plt.ylabel(f'{len(plot_df.index)} Source Neurons', fontsize=16, fontweight='bold')
    
    plt.tight_layout()
    plt.show()
    
    # Return both the dataframe and figure
    return plot_df, fig


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