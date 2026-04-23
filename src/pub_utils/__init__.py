# Import specific functions/classes to make them available at the top level
from .plot import plot_connectome_matrix, plot_similarity_matrix, plot_reciprocal_network, plot_neuron_features
from .core import NeuronFeatures, NeuronInteraction
from .io import handle_pickle, get_file_for_pair, standardize_dataframe, compare_connectomes
from .constants import AllHermNeurons, AllMaleNeurons, SexSharedNeurons, HermSpecificNeurons, MaleSpecificNeurons, AllHermNeuronBlocks
from .assemble import (
    # NT assembly
    get_release_vector,
    get_receptor_matrix,
    assemble_nt_connectome,
    # NPP assembly
    get_npp_release_vector,
    get_npp_receptor_matrix,
    assemble_npp_connectome,
)
from .constrain import (
    load_structural_connectome,
    apply_structural_constraint,
    constrain_assembly,
    get_available_structural_datasets,
)
from .similarity import (
    # Profile functions
    get_sensory_profile,
    get_process_profile,
    get_nt_release_profile,
    get_nt_receptor_profile,
    get_npp_release_profile,
    get_npp_receptor_profile,
    # Similarity computation
    compute_pairwise_similarity,
    compute_release_similarity,
    compute_receptor_similarity,
    compute_sensory_similarity,
    compute_process_similarity,
    compute_all_similarities,
    # Utility functions
    find_similar_neurons,
    compare_similarity_dimensions,
)
from .metrics import (
    compute_sparsity,
    compute_modularity,
    compute_hubness,
    compute_recurrence,
    compute_metrics,
    compute_all_metrics,
    compute_path_length_matrix,
)
from .generate import (
    generate_sparse,
    generate_modular,
    generate_hub,
    generate_reciprocal,
    generate_clustered,
)
from .config import (
    AssemblyConfig,
    ReleaseConfig,
    ReceptorConfig,
    ConstraintConfig,
    OutputConfig,
    load_assembly_config,
    validate_config,
    save_config,
    config_to_dict,
    compute_config_hash,
    create_metadata_sidecar,
)

# Define what is exported when someone does 'from pub_utils import *'
__all__ = [
    "plot_neuron_features",
    "plot_connectome_matrix",
    "plot_similarity_matrix",
    "plot_reciprocal_network",
    "handle_pickle",
    "get_file_for_pair",
    "standardize_dataframe",
    "compare_connectomes",
    "NeuronInteraction",
    "NeuronFeatures",
    "AllHermNeurons",
    "AllMaleNeurons",
    "SexSharedNeurons",
    "HermSpecificNeurons",
    "MaleSpecificNeurons",
    "AllHermNeuronBlocks",
    # Assembly functions
    "get_release_vector",
    "get_receptor_matrix",
    "assemble_nt_connectome",
    "get_npp_release_vector",
    "get_npp_receptor_matrix",
    "assemble_npp_connectome",
    # Structural constraint functions
    "load_structural_connectome",
    "apply_structural_constraint",
    "constrain_assembly",
    "get_available_structural_datasets",
    # Config functions
    "AssemblyConfig",
    "ReleaseConfig",
    "ReceptorConfig",
    "ConstraintConfig",
    "OutputConfig",
    "load_assembly_config",
    "validate_config",
    "save_config",
    "config_to_dict",
    "compute_config_hash",
    "create_metadata_sidecar",
    # Similarity functions
    "get_sensory_profile",
    "get_process_profile",
    "get_nt_release_profile",
    "get_nt_receptor_profile",
    "get_npp_release_profile",
    "get_npp_receptor_profile",
    "compute_pairwise_similarity",
    "compute_release_similarity",
    "compute_receptor_similarity",
    "compute_sensory_similarity",
    "compute_process_similarity",
    "compute_all_similarities",
    "find_similar_neurons",
    "compare_similarity_dimensions",
    # Network metrics
    "compute_sparsity",
    "compute_modularity",
    "compute_hubness",
    "compute_recurrence",
    "compute_metrics",
    "compute_all_metrics",
    "compute_path_length_matrix",
    # Graph generators
    "generate_sparse",
    "generate_modular",
    "generate_hub",
    "generate_reciprocal",
    "generate_clustered",
]

__version__ = "0.1.0"