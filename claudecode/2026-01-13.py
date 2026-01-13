"""
ClaudeCode Session: 2026-01-13

================================================================================
FILES CREATED
================================================================================

data/Bentley2016/
    neuropeptide_receptor_metabotropic_literature.csv  # 13 receptors from Bentley2016
    NT_release_literature.csv                          # 6 markers: cat-2, dat-1, mod-5, tbh-1, tdc-1, tph-1
    NPP_release_literature.csv                         # 31 neuropeptides from Bentley2016

data/RipollSanchez2023/
    NPP_release_sequencing.csv                         # 108 neuropeptides, 302 neurons (109 columns)
    NPP_receptor_all_sequencing.csv                    # 138 GPCRs, 296 neurons with expression (139 columns)

data/Bentley2016/
    NT_receptor_all_literature.csv                     # 17 receptors (merged dopamine, serotonin, octopamine, tyramine)

================================================================================
"""

import pandas as pd
import re
from pathlib import Path

DATA_DIR = Path(__file__).parent.parent / "data"


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


def create_neuropeptide_receptor_metabotropic_file() -> None:
    """
    Convert Bentley2016 neuropeptide_receptor_expression.csv from long format
    to wide format (Dag2023 style).

    Input:
        data/Bentley2016/neuropeptide_receptor_expression.csv (long format)
            Columns: receptor, neuron, wormbase_id, reference

    Output:
        data/Bentley2016/neuropeptide_receptor_metabotropic_literature.csv
            13 receptors: ckr-2, egl-6, frpr-4, npr-1, npr-11, npr-17, npr-2,
                          npr-3, npr-4, npr-5, ntr-1, ntr-2, pdfr-1

    Processing:
        - Standardizes neuronID format (AS1 -> AS01)
        - Uses canonical 302 hermaphrodite neuron order from neuroanatomy.csv
        - Values: 1 = expressed, 0 = not expressed
    """
    # Load the source data
    npp_df = pd.read_csv(DATA_DIR / "Bentley2016/neuropeptide_receptor_expression.csv")

    # Load neuronID order from neuroanatomy.csv
    neuroanatomy = pd.read_csv(DATA_DIR / "RipollSanchez2023/neuroanatomy.csv")
    neuron_order = neuroanatomy['neuronID'].tolist()

    # Get unique receptors
    receptors = sorted(npp_df['receptor'].unique())

    # Standardize neuron names
    npp_df['neuronID'] = npp_df['neuron'].apply(standardize_neuron_id)

    # Create expression column (1 for expressed)
    npp_df['expressed'] = 1

    # Pivot to wide format
    pivot_df = npp_df.pivot_table(
        index='neuronID',
        columns='receptor',
        values='expressed',
        aggfunc='max',
        fill_value=0
    ).reset_index()

    # Reorder columns: neuronID first, then receptors in sorted order
    cols = ['neuronID'] + receptors
    pivot_df = pivot_df[cols]

    # Create full dataframe with all neurons
    full_df = pd.DataFrame({'neuronID': neuron_order})
    full_df = full_df.merge(pivot_df, on='neuronID', how='left')

    # Fill NaN with 0
    full_df = full_df.fillna(0).astype({r: int for r in receptors})

    # Save to CSV
    output_path = DATA_DIR / "Bentley2016/neuropeptide_receptor_metabotropic_literature.csv"
    full_df.to_csv(output_path, index=False)

    # Print summary
    expressed_neurons = (full_df[receptors].sum(axis=1) > 0).sum()
    print(f"Neuropeptide receptors: {len(receptors)} receptors, {expressed_neurons} neurons with expression")
    print(f"  Receptors: {receptors}")
    print(f"  Saved to: {output_path}")


def create_nt_release_literature_file() -> None:
    """
    Convert Bentley2016 monoamine_expression.csv from long format
    to wide format (Wang2024 NT_release_reporter style).

    Input:
        data/Bentley2016/monoamine_expression.csv (long format)
            Columns: monoamine, marker, neuron, wormbase_ID, excluded, excluded_reason, reference

    Output:
        data/Bentley2016/NT_release_literature.csv
            6 markers: cat-2, dat-1, mod-5, tbh-1, tdc-1, tph-1

    Processing:
        - Filters out excluded=1 entries (8 rows with weak/conditional expression)
        - Standardizes neuronID format (AS1 -> AS01)
        - Uses canonical 302 hermaphrodite neuron order from neuroanatomy.csv
        - Values: 1 = expressed, 0 = not expressed
    """
    # Load the source data
    ma_df = pd.read_csv(DATA_DIR / "Bentley2016/monoamine_expression.csv")

    # Load neuronID order from neuroanatomy.csv
    neuroanatomy = pd.read_csv(DATA_DIR / "RipollSanchez2023/neuroanatomy.csv")
    neuron_order = neuroanatomy['neuronID'].tolist()

    # Get unique markers
    markers = sorted(ma_df['marker'].unique())

    # Filter out excluded entries
    df = ma_df[ma_df['excluded'] == 0].copy()

    # Standardize neuron names
    df['neuronID'] = df['neuron'].apply(standardize_neuron_id)

    # Create expression column (1 for expressed)
    df['expressed'] = 1

    # Pivot to wide format
    pivot_df = df.pivot_table(
        index='neuronID',
        columns='marker',
        values='expressed',
        aggfunc='max',
        fill_value=0
    ).reset_index()

    # Reorder columns: neuronID first, then markers in sorted order
    cols = ['neuronID'] + markers
    pivot_df = pivot_df[cols]

    # Create full dataframe with all neurons
    full_df = pd.DataFrame({'neuronID': neuron_order})
    full_df = full_df.merge(pivot_df, on='neuronID', how='left')

    # Fill NaN with 0
    full_df = full_df.fillna(0).astype({m: int for m in markers})

    # Save to CSV
    output_path = DATA_DIR / "Bentley2016/NT_release_literature.csv"
    full_df.to_csv(output_path, index=False)

    # Print summary
    expressed_neurons = (full_df[markers].sum(axis=1) > 0).sum()
    print(f"NT release markers: {len(markers)} markers, {expressed_neurons} neurons with expression")
    print(f"  Markers: {markers}")
    print(f"  Saved to: {output_path}")


def create_npp_release_literature_file() -> None:
    """
    Convert Bentley2016 neuropeptide_expression.csv from long format
    to wide format (similar to NT_release_reporter style).

    Input:
        data/Bentley2016/neuropeptide_expression.csv (long format)
            Columns: neuropeptide, neuron, wormbase_id, reference

    Output:
        data/Bentley2016/NPP_release_literature.csv
            31 neuropeptides: flp-1, flp-2, flp-4, flp-5, flp-7, flp-10, flp-11,
                              flp-13, flp-15, flp-17, flp-18, flp-21, flp-22,
                              nlp-1, nlp-3, nlp-12, nlp-13, nlp-24, ntc-1,
                              pdf-1 through pdf-12

    Processing:
        - Standardizes neuronID format (AS1 -> AS01, but I3 stays I3)
        - Uses canonical 302 hermaphrodite neuron order from neuroanatomy.csv
        - Values: 1 = expressed, 0 = not expressed
    """
    # Load the source data
    npp_df = pd.read_csv(DATA_DIR / "Bentley2016/neuropeptide_expression.csv")

    # Load neuronID order from neuroanatomy.csv
    neuroanatomy = pd.read_csv(DATA_DIR / "RipollSanchez2023/neuroanatomy.csv")
    neuron_order = neuroanatomy['neuronID'].tolist()

    # Get unique neuropeptides
    neuropeptides = sorted(npp_df['neuropeptide'].unique())

    # Standardize neuron names
    npp_df['neuronID'] = npp_df['neuron'].apply(standardize_neuron_id)

    # Create expression column (1 for expressed)
    npp_df['expressed'] = 1

    # Pivot to wide format
    pivot_df = npp_df.pivot_table(
        index='neuronID',
        columns='neuropeptide',
        values='expressed',
        aggfunc='max',
        fill_value=0
    ).reset_index()

    # Reorder columns: neuronID first, then neuropeptides in sorted order
    cols = ['neuronID'] + neuropeptides
    pivot_df = pivot_df[cols]

    # Create full dataframe with all neurons
    full_df = pd.DataFrame({'neuronID': neuron_order})
    full_df = full_df.merge(pivot_df, on='neuronID', how='left')

    # Fill NaN with 0
    full_df = full_df.fillna(0).astype({n: int for n in neuropeptides})

    # Save to CSV
    output_path = DATA_DIR / "Bentley2016/NPP_release_literature.csv"
    full_df.to_csv(output_path, index=False)

    # Print summary
    expressed_neurons = (full_df[neuropeptides].sum(axis=1) > 0).sum()
    print(f"NPP release: {len(neuropeptides)} neuropeptides, {expressed_neurons} neurons with expression")
    print(f"  Saved to: {output_path}")


def convert_ripollsanchez_per_neuron_to_wide(
    input_filename: str,
    output_filename: str,
    entity_type: str = "neuropeptide"
) -> None:
    """
    Convert RipollSanchez2023 per-neuron expression files from column-per-neuron format
    to standard wide format.

    Input format (NPP_per_neuron.csv, GPCR_per_neuron.csv):
        - Row 1: neuronIDs as column headers (302 neurons)
        - Rows 2+: Each cell contains entity name if expressed at that rank, empty if not
        - Multiple rows per neuron indicate ranked expression

    Output format:
        - Column 1: neuronID
        - Remaining columns: one per entity (neuropeptide or GPCR)
        - Values: 1 = expressed, 0 = not expressed

    Args:
        input_filename: Name of input file in RipollSanchez2023 directory
        output_filename: Name of output file in RipollSanchez2023 directory
        entity_type: Description for logging ("neuropeptide" or "GPCR")
    """
    # Load neuronID order from neuroanatomy.csv
    neuroanatomy = pd.read_csv(DATA_DIR / "RipollSanchez2023/neuroanatomy.csv")
    neuron_order = neuroanatomy['neuronID'].tolist()

    # Load per-neuron file (column = neuronID, rows = ranked expression)
    input_df = pd.read_csv(DATA_DIR / f"RipollSanchez2023/{input_filename}")

    # Collect all unique entities and which neurons express them
    all_entities = set()
    neuron_to_entities = {}
    for col in input_df.columns:
        neuron_to_entities[col] = set()
        for val in input_df[col].dropna():
            if str(val).strip():
                all_entities.add(str(val).strip())
                neuron_to_entities[col].add(str(val).strip())

    all_entities = sorted(all_entities)

    # Create output dataframe
    output_data = {'neuronID': neuron_order}
    for entity in all_entities:
        output_data[entity] = []
        for neuron in neuron_order:
            if neuron in neuron_to_entities and entity in neuron_to_entities[neuron]:
                output_data[entity].append(1)
            else:
                output_data[entity].append(0)

    output_df = pd.DataFrame(output_data)
    output_path = DATA_DIR / f"RipollSanchez2023/{output_filename}"
    output_df.to_csv(output_path, index=False)

    # Summary
    expressed_neurons = (output_df[all_entities].sum(axis=1) > 0).sum()
    print(f"{output_filename}:")
    print(f"  {len(all_entities)} {entity_type}s, {expressed_neurons} neurons with expression")
    print(f"  Columns: {output_df.shape[1]} (1 neuronID + {len(all_entities)} {entity_type}s)")
    print(f"  Saved to: {output_path}")


def create_nt_receptor_all_literature_file() -> None:
    """
    Horizontally concatenate monoamine receptor files into a single table.

    Input files:
        data/Bentley2016/dopamine_receptor_all_reporter.csv (5 receptors)
        data/Bentley2016/serotonin_receptor_all_literature.csv (5 receptors)
        data/Bentley2016/octopamine_receptor_all_literature.csv (3 receptors)
        data/Bentley2016/tyramine_receptor_all_literature.csv (4 receptors)

    Output:
        data/Bentley2016/NT_receptor_all_literature.csv
            17 receptors: dop-1, dop-2, dop-3, dop-4, lgc-53, lgc-55, mod-1,
                          octr-1, ser-1, ser-2, ser-3, ser-4, ser-5, ser-6,
                          ser-7, tyra-2, tyra-3
    """
    bentley_dir = DATA_DIR / "Bentley2016"

    # Load all four receptor files
    dopamine_df = pd.read_csv(bentley_dir / 'dopamine_receptor_all_reporter.csv')
    serotonin_df = pd.read_csv(bentley_dir / 'serotonin_receptor_all_literature.csv')
    octopamine_df = pd.read_csv(bentley_dir / 'octopamine_receptor_all_literature.csv')
    tyramine_df = pd.read_csv(bentley_dir / 'tyramine_receptor_all_literature.csv')

    # Merge on neuronID
    merged_df = dopamine_df.merge(serotonin_df, on='neuronID', how='outer')
    merged_df = merged_df.merge(octopamine_df, on='neuronID', how='outer')
    merged_df = merged_df.merge(tyramine_df, on='neuronID', how='outer')

    # Sort columns: neuronID first, then receptors alphabetically
    receptor_cols = sorted([col for col in merged_df.columns if col != 'neuronID'])
    merged_df = merged_df[['neuronID'] + receptor_cols]

    # Save
    output_path = bentley_dir / 'NT_receptor_all_literature.csv'
    merged_df.to_csv(output_path, index=False)

    print(f"NT_receptor_all_literature.csv:")
    print(f"  {len(receptor_cols)} receptors, {merged_df.shape[0]} neurons")
    print(f"  Receptors: {receptor_cols}")
    print(f"  Saved to: {output_path}")


if __name__ == "__main__":
    create_neuropeptide_receptor_metabotropic_file()
    create_nt_release_literature_file()
    create_npp_release_literature_file()
    create_nt_receptor_all_literature_file()

    # RipollSanchez2023 conversions
    convert_ripollsanchez_per_neuron_to_wide(
        "NPP_per_neuron.csv",
        "NPP_release_sequencing.csv",
        "neuropeptide"
    )
    convert_ripollsanchez_per_neuron_to_wide(
        "GPCR_per_neuron.csv",
        "NPP_receptor_all_sequencing.csv",
        "GPCR"
    )
