"""
ClaudeCode Session: 2026-01-12
Commit: 21e65a6 - massive refactoring of raw data for Bentley2016, Fenyves2020 and Wang2024

================================================================================
CURATED DATA FILES CREATED WITH CLAUDE ASSISTANCE
================================================================================

data/Fenyves2020/
    acetylcholine_receptor_ionotropic_sequencing.csv
    gaba_receptor_ionotropic_sequencing.csv
    glutamate_receptor_ionotropic_sequencing.csv

data/HobertLab/
    acetylcholine_receptor_metabotropic_reporter.csv   # neuronID x gar-1/gar-2/gar-3
    gaba_receptor_all_reporter.csv                     # neuronID x GABA receptor genes

data/Wang2024/
    NT_release_reporter.csv          # hermaphrodite NT release (reporter-based)
    NT_release_reporter_male.csv     # male NT release (reporter-based)
    NT_release_staining.csv          # hermaphrodite NT release (staining-based)
    NT_release_staining_male.csv     # male NT release (staining-based)

data/Bentley2016/
    *.csv                               # TSV to CSV conversions
    serotonin_receptor_all_literature.csv   # 5 receptors: mod-1, ser-1, ser-4, ser-5, ser-7
    dopamine_receptor_all_reporter.csv      # 5 receptors: dop-1, dop-2, dop-3, dop-4, lgc-53
    octopamine_receptor_all_literature.csv  # 3 receptors: octr-1, ser-3, ser-6
    tyramine_receptor_all_literature.csv    # 4 receptors: lgc-55, ser-2, tyra-2, tyra-3

================================================================================
CODE CHANGES: src/pub_utils/io.py
================================================================================
"""

import pickle
import pandas as pd
from pathlib import Path

# Default path to neuroanatomy mapping file
# __file__ is src/pub_utils/io.py, so .parent.parent.parent gets to project root
_DEFAULT_MAPPING_PATH = Path(__file__).parent.parent.parent / 'data' / 'RipollSanchez2023' / 'neuroanatomy.csv'


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


# ==============================================================================
# SESSION 2: Bentley2016 Monoamine Receptor Expression Files
# ==============================================================================
#
# Files created:
#   data/Bentley2016/*.csv (TSV to CSV conversions)
#   data/Bentley2016/serotonin_receptor_all_literature.csv
#   data/Bentley2016/dopamine_receptor_all_reporter.csv
#   data/Bentley2016/octopamine_receptor_all_literature.csv
#   data/Bentley2016/tyramine_receptor_all_literature.csv
# ==============================================================================

import re
import subprocess

DATA_DIR = Path(__file__).parent.parent / "data"


def convert_tsv_to_csv(directory: Path) -> None:
    """Convert all TSV files in a directory to CSV format."""
    for tsv_file in directory.glob("*.tsv"):
        csv_file = tsv_file.with_suffix(".csv")
        subprocess.run(
            f"sed 's/\\t/,/g' \"{tsv_file}\" > \"{csv_file}\"",
            shell=True,
            check=True
        )
        print(f"Converted {tsv_file.name} -> {csv_file.name}")


def standardize_neuron_id(name: str) -> str:
    """
    Standardize neuronID format (AS1 -> AS01, DD1 -> DD01, etc.).

    Only zero-pads neurons with 2+ letter prefixes (AS, DA, DD, VA, VB, VD, VC, etc.)
    Single-letter prefixes (I, M) are not zero-padded to match canonical names.
    """
    match = re.match(r'^([A-Z]{2,})(\d+)$', name)  # Only 2+ letter prefixes
    if match:
        prefix, num = match.groups()
        return f"{prefix}{int(num):02d}"
    return name


def create_monoamine_receptor_expression_files() -> None:
    """
    Convert Bentley2016 monoamine_receptor_expression.csv from long format
    to wide format (Dag2023 style), creating one file per monoamine.

    Input:
        data/Bentley2016/monoamine_receptor_expression.csv (long format)
            Columns: monoamine, receptor, neuron, wormbase_id, excluded, excluded_reason, reference

    Output files:
        data/Bentley2016/serotonin_receptor_all_literature.csv (5 receptors: mod-1, ser-1, ser-4, ser-5, ser-7)
        data/Bentley2016/dopamine_receptor_all_reporter.csv (5 receptors: dop-1, dop-2, dop-3, dop-4, lgc-53)
        data/Bentley2016/octopamine_receptor_all_literature.csv (3 receptors: octr-1, ser-3, ser-6)
        data/Bentley2016/tyramine_receptor_all_literature.csv (4 receptors: lgc-55, ser-2, tyra-2, tyra-3)

    Processing:
        - Filters out excluded=1 entries (67 rows with hypothetical/weak expression)
        - Standardizes neuronID format (AS1 -> AS01)
        - Uses canonical 302 hermaphrodite neuron order from neuroanatomy.csv
        - Values: 1 = expressed, 0 = not expressed
    """
    # Load the source data
    bentley_df = pd.read_csv(DATA_DIR / "Bentley2016/monoamine_receptor_expression.csv")

    # Load neuronID order from neuroanatomy.csv
    neuroanatomy = pd.read_csv(DATA_DIR / "RipollSanchez2023/neuroanatomy.csv")
    neuron_order = neuroanatomy['neuronID'].tolist()

    # Filter out excluded entries
    df = bentley_df[bentley_df['excluded'] == 0].copy()

    # Standardize neuron names
    df['neuronID'] = df['neuron'].apply(standardize_neuron_id)

    # Process each monoamine with output filename mapping
    monoamine_files = {
        'serotonin': 'serotonin_receptor_all_literature.csv',
        'dopamine': 'dopamine_receptor_all_reporter.csv',
        'octopamine': 'octopamine_receptor_all_literature.csv',
        'tyramine': 'tyramine_receptor_all_literature.csv',
    }

    for monoamine in monoamine_files:
        # Filter for this monoamine
        ma_df = df[df['monoamine'] == monoamine].copy()

        # Get unique receptors for this monoamine
        receptors = sorted(ma_df['receptor'].unique())

        # Create expression column (1 for expressed)
        ma_df['expressed'] = 1

        # Pivot to wide format
        pivot_df = ma_df.pivot_table(
            index='neuronID',
            columns='receptor',
            values='expressed',
            aggfunc='max',  # In case of duplicates, take max
            fill_value=0
        ).reset_index()

        # Ensure all receptors are present as columns
        for receptor in receptors:
            if receptor not in pivot_df.columns:
                pivot_df[receptor] = 0

        # Reorder columns: neuronID first, then receptors in sorted order
        cols = ['neuronID'] + receptors
        pivot_df = pivot_df[cols]

        # Create full dataframe with all neurons
        full_df = pd.DataFrame({'neuronID': neuron_order})
        full_df = full_df.merge(pivot_df, on='neuronID', how='left')

        # Fill NaN with 0
        full_df = full_df.fillna(0).astype({r: int for r in receptors})

        # Save to CSV
        output_path = DATA_DIR / f"Bentley2016/{monoamine_files[monoamine]}"
        full_df.to_csv(output_path, index=False)

        # Print summary
        expressed_neurons = (full_df[receptors].sum(axis=1) > 0).sum()
        print(f"{monoamine}: {len(receptors)} receptors, {expressed_neurons} neurons with expression")
        print(f"  Receptors: {receptors}")
        print(f"  Saved to: {output_path}")
        print()


if __name__ == "__main__":
    # Convert TSV to CSV
    print("Converting TSV files to CSV...")
    convert_tsv_to_csv(DATA_DIR / "Bentley2016")
    print()

    # Create monoamine receptor expression files
    print("Creating monoamine receptor expression files...")
    create_monoamine_receptor_expression_files()
    print("Done!")
