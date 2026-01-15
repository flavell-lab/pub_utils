# Reproducibility Test: NPP_RipollSanchez2023

**Date**: 2026-01-15
**Tested by**: Claude Code (assemble.py validation)

## Objective

Validate that `assemble_npp_connectome()` produces results consistent with the preassembled RipollSanchez2023 neuropeptide connectomes.

## Data Sources

### Preassembled (reference)
- **Files**: `connectomes/preassembled/NPP_RipollSanchez2023_{shortRange,midRange,longRange}.csv`
- **Origin**: OpenWorm Connectome Toolbox (`cect.RipollSanchezShortRangeReader`, etc.)

| Range | Sum | Non-zero entries |
|-------|-----|------------------|
| shortRange | 83,634 | 31,417 |
| midRange | 107,429 | 40,425 |
| longRange | 145,834 | 53,558 |

### Assembly inputs
- **Release data**: `data/RipollSanchez2023/NPP_release_sequencing.csv` (108 neuropeptides)
- **Receptor data**: `data/RipollSanchez2023/NPP_receptor_all_sequencing.csv` (138 GPCRs)
- **Pairing info**: `data/RipollSanchez2023/NPP_receptor_info.csv` (92 pairs, 49 unique NPPs)

## Test: Aggregate all NPP connectomes

Assembled all 49 NPPs from RipollSanchez2023 pairing info and summed the count matrices.

```python
aggregate = sum([
    assemble_npp_connectome(npp, output_format='count', ...)
    for npp in all_npps
])
```

### Result

| Metric | My Assembly | longRange | midRange | shortRange |
|--------|-------------|-----------|----------|------------|
| Sum | **145,834** | **145,834** | 107,429 | 83,634 |
| Non-zero | **53,558** | **53,558** | 40,425 | 31,417 |

## Conclusion

**Assembly logic validated.** My aggregate exactly matches the **longRange** preassembled connectome.

This is expected because:
1. My assembly includes ALL NPP-receptor pairs without distance constraints
2. The "long-range" model in RipollSanchez2023 also has no distance constraints
3. The "short-range" and "mid-range" models apply physical proximity filters based on neuronal anatomy

## Notes on range models

The RipollSanchez2023 paper models neuropeptide signaling at three ranges:
- **Short-range**: Neuropeptides only reach nearby neurons (anatomical proximity)
- **Mid-range**: Intermediate diffusion distance
- **Long-range**: Neuropeptides can reach any receptor-expressing neuron (no distance constraint)

Our assembly function currently implements the long-range model. Adding distance constraints would require anatomical proximity data not currently in our data assets.

## Assembly details

- **49 NPPs assembled**: flp-1 through flp-28, nlp-* series, ins-* series, pdf-*, etc.
- **0 NPPs skipped**: All NPPs in pairing info had matching release and receptor data
