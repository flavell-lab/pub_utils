"""
Claude Code Session Log - 2026-01-15
=====================================

Task: Recreate serotonin_receptor_all_reporter.csv from Dag2023 source data

SOURCE FILE
-----------
data/Dag2023/5htr_expression_dv_final.csv
- Format: Neurons (abbreviated class names), ser1, ser5, ser7, lgc50, ser4, mod1
- Values: 1 for positive expression, empty for no data
- Contains 64 neuron entries with abbreviated names (e.g., RMD, RMDD, CEPD)

TARGET FILE
-----------
data/Dag2023/serotonin_receptor_all_reporter.csv
- Format: neuronID, ser-1, ser-5, ser-7, lgc-50, ser-4, mod-1
- Values: 1 for positive, 0 for negative/unstudied
- Must cover all 302 hermaphrodite neurons with proper L/R suffixes

TRANSFORMATION LOGIC
--------------------
1. Column renaming: ser1 → ser-1, ser5 → ser-5, etc. (add hyphens)
2. Neuron name expansion: Use neuroanatomy.csv to map neuronClass → neuronID list
   - Primary mapping: neuroanatomy.csv neuronClass → list of neuronIDs
   - Subclass mapping: Automatically derived for D/V subclasses not in neuroanatomy
     (e.g., RMDD → RMDDL/RMDDR, CEPV → CEPVL/CEPVR)
   - Motor neuron format: DD1 → DD01 (zero-padded format)
3. Fill missing neurons with 0 (unstudied)
4. Replace empty values with 0

OUTPUT STATISTICS
-----------------
- 109 neurons with any receptor expression
- ser-1: 45 neurons
- ser-5: 35 neurons
- ser-7: 20 neurons
- lgc-50: 34 neurons
- ser-4: 56 neurons
- mod-1: 36 neurons

EXECUTION CODE
--------------
"""

import pandas as pd
import numpy as np
import sys
sys.path.insert(0, '/storage/fs/store1/shared/pub_utils/src')
from pub_utils.constants import AllHermNeurons

# Load neuroanatomy to build class → neuronIDs mapping
neuroanatomy_path = '/storage/fs/store1/shared/pub_utils/data/RipollSanchez2023/neuroanatomy.csv'
neuroanatomy = pd.read_csv(neuroanatomy_path)

# Build mapping from neuronClass to list of neuronIDs
class_to_neurons = neuroanatomy.groupby('neuronClass')['neuronID'].apply(list).to_dict()

# Add subclass mappings for neurons that use D/V suffixes in source data
# These are more specific than the broad class names in neuroanatomy
# Pattern: sourceClass + D/V → neuronIDs containing D/V
subclass_expansions = {}
for neuron_id in neuroanatomy['neuronID']:
    # Check for D/V subclass patterns (e.g., RMDDL → subclass RMDD)
    if len(neuron_id) >= 4:
        # Try extracting subclass (e.g., RMDDL → RMDD, CEPVR → CEPV)
        for suffix_len in [2, 1]:  # Try removing L/R first, then just last char
            if neuron_id.endswith('L') or neuron_id.endswith('R'):
                potential_subclass = neuron_id[:-1]
                if potential_subclass not in class_to_neurons:
                    if potential_subclass not in subclass_expansions:
                        subclass_expansions[potential_subclass] = []
                    if neuron_id not in subclass_expansions[potential_subclass]:
                        subclass_expansions[potential_subclass].append(neuron_id)

# Also add numeric motor neuron patterns (DD1 → DD01)
subclass_expansions['DD1'] = ['DD01']
subclass_expansions['VB02'] = ['VB02']  # Already correct format

# Merge subclass expansions into class_to_neurons
class_to_neurons.update(subclass_expansions)

# Read source file
source_path = '/storage/fs/store1/shared/pub_utils/data/Dag2023/5htr_expression_dv_final.csv'
df = pd.read_csv(source_path)

# Column renaming: add hyphens to gene names
column_map = {
    'Neurons': 'neuronClass',  # Rename to neuronClass since these are class names
    'ser1': 'ser-1',
    'ser5': 'ser-5',
    'ser7': 'ser-7',
    'lgc50': 'lgc-50',
    'ser4': 'ser-4',
    'mod1': 'mod-1'
}
df = df.rename(columns=column_map)

# Replace empty/NaN with 0
df = df.fillna(0)

# Convert to int where possible
receptors = ['ser-1', 'ser-5', 'ser-7', 'lgc-50', 'ser-4', 'mod-1']
for col in receptors:
    df[col] = df[col].astype(int)

# Build expanded dataframe using neuroanatomy mapping
expanded_rows = []

for _, row in df.iterrows():
    neuron_class = row['neuronClass']
    values = {r: row[r] for r in receptors}

    if neuron_class in class_to_neurons:
        # Expand to all neuronIDs in this class
        for neuron_id in class_to_neurons[neuron_class]:
            expanded_rows.append({'neuronID': neuron_id, **values})
    else:
        # Try direct match (in case source uses neuronID instead of class)
        if neuron_class in AllHermNeurons:
            expanded_rows.append({'neuronID': neuron_class, **values})
        else:
            print(f"Warning: '{neuron_class}' not found in neuroanatomy class_to_neurons or AllHermNeurons")

expanded_df = pd.DataFrame(expanded_rows)

# Create final dataframe with all 302 neurons, default 0
final_df = pd.DataFrame(0, index=AllHermNeurons, columns=receptors)
final_df.index.name = 'neuronID'

# Fill in values from expanded data
for _, row in expanded_df.iterrows():
    neuron = row['neuronID']
    if neuron in final_df.index:
        for r in receptors:
            final_df.loc[neuron, r] = row[r]
    else:
        print(f"Warning: {neuron} not in AllHermNeurons")

# Reset index to make neuronID a column
final_df = final_df.reset_index()

# Save
output_path = '/storage/fs/store1/shared/pub_utils/data/Dag2023/serotonin_receptor_all_reporter.csv'
final_df.to_csv(output_path, index=False)

print(f"Created {output_path}")
print(f"Shape: {final_df.shape}")
print(f"Neurons with any receptor: {(final_df[receptors].sum(axis=1) > 0).sum()}")
print(f"Receptor coverage:")
for r in receptors:
    print(f"  {r}: {(final_df[r] > 0).sum()} neurons")

"""
CHANGES TO data/assets.json
---------------------------
Fixed serotonin receptor entry - changed method from "literature" to "reporter":

Before:
  "serotonin": {
    "all": {
      "literature": {
        "Dag2023": "data/Dag2023/serotonin_receptor_all_reporter.csv"
      }
    }
  }

After:
  "serotonin": {
    "all": {
      "reporter": {
        "Dag2023": "data/Dag2023/serotonin_receptor_all_reporter.csv"
      }
    }
  }


CHANGES TO notebook/assemble_connectomes.ipynb
----------------------------------------------
Fixed cell dcicmr9qt5p - changed receptor_sources from 'literature:Dag2023' to 'reporter:Dag2023'
and updated description text from "literature-based" to "reporter-based".


BUG FIX TO src/pub_utils/plot.py
--------------------------------
Changed line 20 from:
    actual_max = int(plot_df.max().max())
To:
    actual_max = int(np.nanmax(plot_df.values))

This prevents ValueError when plot_df contains NaN values, as np.nanmax ignores NaNs.


SESSION NOTES
-------------
- Root cause of serotonin all-NaN issue: assets.json had "literature" as method key but
  notebook was calling receptor_sources=['reporter:Dag2023']. The code couldn't find the path.
- The file path in assets.json was already correct (serotonin_receptor_all_reporter.csv),
  only the method key was wrong.
- release_markers=['synthesis'] for serotonin resolves to 'tph-1' gene, which exists in
  data/Bentley2016/NT_release_literature.csv - so release side was working correctly.
- Used neuroanatomy.csv to derive neuron class → neuronID mappings instead of hardcoding
  expansion rules. This is more robust and maintainable.


CHANGES TO src/pub_utils/io.py
------------------------------
Added compare_connectomes() function for comparing two connectome CSVs:

    result = pu.compare_connectomes(path_a, path_b, names=None, threshold=0)

Returns dict with:
    - 'summary': DataFrame with metrics (connections count, overlap, Jaccard, correlation)
    - 'overlap': DataFrame of connections in both
    - 'only_{name_a}': DataFrame of connections only in A
    - 'only_{name_b}': DataFrame of connections only in B
    - 'diff': DataFrame of A - B values
    - 'correlation': Pearson correlation value

Example output:
                    Metric     Value
    Neurons (common)          302
    Connections in A         3709
    Connections in B         2194
    Overlap (in both)        2002
    Only in A                1707
    Only in B                 192
    Jaccard similarity      0.5132
    Pearson correlation     0.6232


CHANGES TO src/pub_utils/__init__.py
------------------------------------
Added compare_connectomes to imports and __all__ exports.


CHANGES TO notebook/assemble_connectomes.ipynb
----------------------------------------------
Updated NPP example (Example 3) to follow same pattern as NT examples:
- Full assembly with metadata saving
- Binary and count format outputs
- Per-receptor connectome outputs
- Removed outdated flp1_conn and flp1_per_receptor examples

Fixed section and example numbering:
- Sections: 1 (NT), 2 (NPP), 3 (Lower-Level), 4 (Aggregating)
- Examples: 1 (Dopamine), 2 (Serotonin), 3 (FLP-1), 4 (Inspect), 5 (Aggregate)
"""
