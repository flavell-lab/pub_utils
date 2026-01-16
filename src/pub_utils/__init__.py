# Import specific functions/classes to make them available at the top level
from .plot import plot_connectome_matrix, plot_reciprocal_network, plot_neuron_features
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

# Define what is exported when someone does 'from pub_utils import *'
__all__ = [
    "plot_neuron_features",
    "plot_connectome_matrix",
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
]

__version__ = "0.1.0"