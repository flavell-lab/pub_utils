# Reproducibility Test: dopamine_Bentley2016

**Date**: 2026-01-15
**Tested by**: Claude Code (assemble.py validation)

## Objective

Validate that `assemble_nt_connectome()` produces results consistent with the preassembled `dopamine_Bentley2016.csv` connectome from the OpenWorm `cect` library.

## Data Sources

### Preassembled (reference)
- **File**: `connectomes/preassembled/molecular/dopamine_Bentley2016.csv`
- **Origin**: OpenWorm Connectome Toolbox (`cect.WormNeuroAtlasMAReader`)
- **Connections**: 1160

### Assembly inputs
- **Release data**: `data/Bentley2016/NT_release_literature.csv` (cat-2 marker)
- **Receptor data (literature)**: `data/Bentley2016/NT_receptor_all_literature.csv`
- **Receptor data (reporter)**: `data/Bentley2016/dopamine_receptor_all_reporter.csv`

## Test: Assembly using Bentley2016 receptor definitions

Used Bentley2016's own receptor columns (dop-1, dop-2, dop-3, dop-4, lgc-53) from combined literature + reporter data with OR gate.

### Result
| Metric | My Assembly | Preassembled | Difference |
|--------|-------------|--------------|------------|
| Total connections | 1176 | 1160 | +16 |
| Source neurons | 8 | 8 | 0 |
| Targets per source | 147 | 145 | +2 |

### Root cause
The +16 connections (2 extra targets × 8 sources) are to **CANL** and **CANR**.

The `cect` library explicitly excludes CAN neurons:
> "*This version of the NeuroAtlas does not include the CAN neurons.*"

## Conclusion

**Assembly logic validated.** When using matching input data, the only discrepancy is CAN neuron exclusion by the `cect` library.

## Receptor expression summary

| Receptor | Literature file | Reporter file | Notes |
|----------|-----------------|---------------|-------|
| dop-1 | 72 neurons | present | GPCR |
| dop-2 | 20 neurons | present | GPCR |
| dop-3 | 83 neurons | present | GPCR |
| dop-4 | 10 neurons | present (more neurons) | GPCR |
| lgc-53 | not present | 18 neurons | Ion channel |
