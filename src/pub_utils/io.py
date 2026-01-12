# import h5py
# import zarr
import pickle
import pandas as pd
from pathlib import Path

# Default path to neuroanatomy mapping file
# __file__ is src/pub_utils/io.py, so .parent.parent.parent gets to project root
_DEFAULT_MAPPING_PATH = Path(__file__).parent.parent.parent / 'data' / 'RipollSanchez2023' / 'neuroanatomy.csv'


def handle_pickle(data=None, filename="data.pkl", mode="save"):
    if mode == "save":
        if data is None:
            raise ValueError("You must provide data to save.")
        with open(filename, "wb") as f:
            pickle.dump(data, f)
        print(f"Successfully saved to {filename}")
        
    elif mode == "load":
        try:
            with open(filename, "rb") as f:
                return pickle.load(f)
        except FileNotFoundError:
            print(f"Error: The file {filename} was not found.")
            return None
    else:
        print("Invalid mode. Please use 'save' or 'load'.")


def create_neuropeptide_mapping():
    """
    Maps neuropeptide-receptor pairs to their corresponding file numbers (001-092).
    Returns a case-insensitive dictionary with the pair as key and the file number as value.
    """
    pairs = [
        "NLP-40 AEX-2",
        "NLP-12 CKR-1",
        "NLP-12 CKR-2",
        "FLP-10 DMSR-1",
        "FLP-11 DMSR-1",
        "FLP-12 DMSR-1",
        "FLP-13 DMSR-1",
        "FLP-16 DMSR-1",
        "FLP-17 DMSR-1",
        "FLP-25 DMSR-1",
        "FLP-28 DMSR-1",
        "FLP-33 DMSR-1",
        "FLP-4 DMSR-1",
        "FLP-9 DMSR-1",
        "NLP-13 DMSR-2",
        "FLP-14 DMSR-3",
        "FLP-1 DMSR-5",
        "FLP-1 DMSR-6",
        "FLP-1 DMSR-7",
        "FLP-11 DMSR-7",
        "FLP-13 DMSR-7",
        "FLP-16 DMSR-7",
        "FLP-23 DMSR-7",
        "FLP-27 DMSR-7",
        "FLP-4 DMSR-7",
        "FLP-5 DMSR-7",
        "FLP-7 DMSR-7",
        "FLP-9 DMSR-7",
        "FLP-12 DMSR-8",
        "FLP-10 EGL-6",
        "FLP-17 EGL-6",
        "FLP-4 EGL-6",
        "FLP-5 EGL-6",
        "FLP-6 EGL-6",
        "FLP-9 EGL-6",
        "FLP-8 FRPR-15",
        "FLP-3 FRPR-16",
        "FLP-2 FRPR-18",
        "FLP-14 FRPR-19",
        "FLP-8 FRPR-19",
        "FLP-20 FRPR-3",
        "FLP-33 FRPR-6",
        "FLP-1 FRPR-7",
        "FLP-10 FRPR-8",
        "FLP-11 FRPR-8",
        "FLP-12 FRPR-8",
        "FLP-13 FRPR-8",
        "FLP-14 FRPR-8",
        "FLP-16 FRPR-8",
        "FLP-17 FRPR-8",
        "FLP-25 FRPR-8",
        "FLP-32 FRPR-8",
        "FLP-4 FRPR-8",
        "FLP-8 FRPR-8",
        "FLP-9 FRPR-8",
        "FLP-19 FRPR-9",
        "NLP-47 GNRR-1",
        "NLP-23 GNRR-3",
        "NLP-2 GNRR-6",
        "NLP-22 GNRR-6",
        "NLP-23 GNRR-6",
        "NLP-44 NMUR-1",
        "NLP-44 NMUR-2",
        "FLP-21 NPR-1",
        "FLP-3 NPR-10",
        "FLP-34 NPR-11",
        "FLP-21 NPR-2",
        "NLP-72 NPR-22",
        "FLP-15 NPR-3",
        "NLP-64 NPR-32",
        "SNET-1 NPR-34",
        "NLP-10 NPR-35",
        "NLP-17 NPR-37",
        "FLP-14 NPR-39",
        "FLP-8 NPR-39",
        "FLP-18 NPR-4",
        "FLP-14 NPR-40",
        "FLP-8 NPR-40",
        "NLP-13 NPR-41",
        "NLP-3 NPR-42",
        "NLP-17 NPR-43",
        "FLP-18 NPR-5",
        "FLP-27 NPR-8",
        "FLP-8 NTR-1",
        "NTC-1 NTR-1",
        "NLP-37 PDFR-1",
        "NLP-49 SEB-3",
        "NLP-42 SPRR-1",
        "NLP-38 SPRR-2",
        "NLP-58 TKR-1",
        "NLP-58 TKR-2",
        "NLP-54 TRHR-1"
    ]
    
    # Create case-insensitive mapping dictionary
    mapping = {}
    for i, pair in enumerate(pairs, start=1):
        # Format file number with zero padding
        file_num = f"{i:03d}"
        # Store with lowercase key for case-insensitive lookup
        mapping[pair.lower()] = file_num
    
    return mapping


def get_file_for_pair(neuropeptide, receptor):
    """
    Get the file number for a specific neuropeptide-receptor pair (case-insensitive).
    
    Args:
        neuropeptide: The neuropeptide name (e.g., "FLP-10", "flp-10", "Flp-10")
        receptor: The receptor name (e.g., "DMSR-1", "dmsr-1", "Dmsr-1")
    
    Returns:
        The file number (e.g., "004") or None if not found
    """
    mapping = create_neuropeptide_mapping()
    # Convert to lowercase for case-insensitive lookup
    key = f"{neuropeptide} {receptor}".lower()
    return mapping.get(key)


def standardize_dataframe(df, neuron_order, mapping_df='default', verbose=True):
    """
    Standardize a square connectome matrix to use canonical neuronID naming.

    Args:
        df: Square DataFrame with neuron identifiers as row/column labels.
            If 'Row' column exists, it will be used as the index.
        neuron_order: List of neuronIDs defining the desired row/column order.
        mapping_df: DataFrame containing neuron name mappings, or 'default' to use
            neuroanatomy.csv, or None to skip mapping.
            Must have 'neuronID' column and can have 'neuronClass',
            'figNeuronClass', and/or 'synonym' columns for mapping.
        verbose: If True, print messages about name mappings and missing neurons.

    Returns:
        Standardized DataFrame with neuronID as row/column labels,
        reindexed to match neuron_order.
    """
    plot_df = df.set_index('Row').rename_axis(None) if 'Row' in df.columns else df.copy()

    # Check that it's square
    assert plot_df.shape[0] == plot_df.shape[1], "Not a square matrix"

    # Check that row and column names match
    assert list(plot_df.index) == list(plot_df.columns), "Row and column labels don't match"

    # Load default mapping if requested
    if mapping_df == 'default':
        if _DEFAULT_MAPPING_PATH.exists():
            mapping_df = pd.read_csv(_DEFAULT_MAPPING_PATH)
        else:
            raise FileNotFoundError(f"Default mapping file not found: {_DEFAULT_MAPPING_PATH}")

    # Map row/column names to neuronID
    if mapping_df is not None:
        name_to_id = _build_neuron_mapping(mapping_df)

        # Track alternative name mappings for verbose output
        alt_mappings = []
        for name in plot_df.index:
            if name in name_to_id and name_to_id[name] != name:
                alt_mappings.append((name, name_to_id[name]))

        if verbose and alt_mappings:
            print(f"Mapped {len(alt_mappings)} alternative names to neuronID:")
            for alt_name, neuron_id in alt_mappings:
                print(f"  {alt_name} -> {neuron_id}")

        # Map index and columns to neuronID
        new_index = [name_to_id.get(name, name) for name in plot_df.index]
        new_columns = [name_to_id.get(name, name) for name in plot_df.columns]

        plot_df.index = new_index
        plot_df.columns = new_columns

    # Track missing neurons before reindex
    if verbose:
        current_neurons = set(plot_df.index)
        missing_neurons = [n for n in neuron_order if n not in current_neurons]
        if missing_neurons:
            print(f"Filling {len(missing_neurons)} missing neuronIDs with NaN:")
            for neuron_id in missing_neurons:
                print(f"  {neuron_id} (not found)")

    std_df = plot_df.reindex(index=neuron_order, columns=neuron_order)

    return std_df


def _build_neuron_mapping(mapping_df):
    """
    Build a dictionary mapping alternative neuron names to canonical neuronID.

    Args:
        mapping_df: DataFrame with 'neuronID' and optional columns
            'neuronClass', 'figNeuronClass', 'synonym'.

    Returns:
        Dict mapping any known name to its canonical neuronID.
    """
    if 'neuronID' not in mapping_df.columns:
        raise ValueError("mapping_df must contain 'neuronID' column")

    name_to_id = {}
    mapping_cols = ['neuronClass', 'figNeuronClass', 'synonym']

    for _, row in mapping_df.iterrows():
        neuron_id = row['neuronID']

        # Map neuronID to itself
        name_to_id[neuron_id] = neuron_id

        # Map alternative names to neuronID
        for col in mapping_cols:
            if col in mapping_df.columns:
                alt_name = row.get(col)
                if pd.notna(alt_name) and alt_name != '':
                    # Only map if not already mapped (first occurrence wins)
                    # This handles cases where multiple neuronIDs share a class
                    if alt_name not in name_to_id:
                        name_to_id[alt_name] = neuron_id

    return name_to_id