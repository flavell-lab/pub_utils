# Reproducibility Test: tyramine_Bentley2016

**Date**: 2026-01-15
**Tested by**: Claude Code (assemble.py validation)

## Objective

Validate that `assemble_nt_connectome()` produces results consistent with the preassembled `tyramine_Bentley2016.csv` connectome from the OpenWorm `cect` library.

## Data Sources

### Preassembled (reference)
- **File**: `connectomes/preassembled/tyramine_Bentley2016.csv`
- **Origin**: OpenWorm Connectome Toolbox (`cect.WormNeuroAtlasMAReader`)
- **Connections**: 224

### Assembly inputs
- **Release data**: `data/Bentley2016/NT_release_literature.csv` (tdc-1 marker)
- **Receptor data**: `data/Bentley2016/tyramine_receptor_all_literature.csv`

## Test: Assembly using Bentley2016 receptor definitions

Used Bentley2016's tyramine receptor columns (lgc-55, ser-2, tyra-2, tyra-3) with OR gate across neurons.

### Result
| Metric | My Assembly | Preassembled | Difference |
|--------|-------------|--------------|------------|
| Total connections | 228 | 224 | +4 |
| Source neurons | 2 | 2 | 0 |
| Targets per source | 114 | 112 | +2 |

### Source neurons
- RIML, RIMR (tdc-1 positive)

### Root cause
The +4 connections (2 extra targets × 2 sources) are to **CANL** and **CANR**.

The `cect` library explicitly excludes CAN neurons:
> "*This version of the NeuroAtlas does not include the CAN neurons.*"

## Conclusion

**Assembly logic validated.** When using matching input data, the only discrepancy is CAN neuron exclusion by the `cect` library.

## Receptor expression summary

| Receptor | Neurons expressing | Notes |
|----------|-------------------|-------|
| lgc-55 | 27 | Ion channel |
| ser-2 | 59 | GPCR |
| tyra-2 | 18 | GPCR |
| tyra-3 | 27 | GPCR |
| **Any** | **114** | Union |
