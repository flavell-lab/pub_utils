# import h5py
# import zarr
import pickle


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


def standardize_dataframe(df, neuron_order):
    
    plot_df = df.set_index('Row').rename_axis(None) if 'Row' in df.columns else df.copy()

    # Check that it's square
    assert plot_df.shape[0] == plot_df.shape[1], "Not a square matrix"

    # Check that row and column names match
    assert list(plot_df.index) == list(plot_df.columns), "Row and column labels don't match"
    
    std_df = plot_df.reindex(index=neuron_order, columns=neuron_order)   # Reindex to include all names in col_order, filling missing ones with NaN

    return std_df