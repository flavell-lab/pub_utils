# Reproducibility Test: FLP-1 / DMSR-7 Connectome (RipollSanchez2023)

**Date**: 2026-01-15
**Status**: PASSED

## Objective

Verify that `assemble_npp_connectome()` reproduces the preassembled FLP-1/DMSR-7 molecular connectome from RipollSanchez2023 data.

## Reference File

- `connectomes/preassembled/molecular/FLP-1_DMSR-7_RipollSanchez2023_longRange.csv`

## Test Procedure

### 1. Assembly Parameters

```python
from pub_utils.assemble import assemble_npp_connectome

connectome = assemble_npp_connectome(
    neuropeptide='flp-1',
    release_sources=['sequencing:RipollSanchez2023'],
    receptor_sources=['sequencing:RipollSanchez2023'],
    pairing_source='RipollSanchez2023',
    receptor_gate='or',
    output_format='per_pair'
)

# Extract DMSR-7 specific connectome
dmsr7_connectome = connectome['dmsr-7']
```

### 2. Data Sources Used

| Data Type | Source Key | File Path |
|-----------|------------|-----------|
| Pairing Info | RipollSanchez2023 | `data/RipollSanchez2023/NPP_receptor_info.csv` |
| Release Data | sequencing:RipollSanchez2023 | `data/RipollSanchez2023/NPP_release_sequencing.csv` |
| Receptor Data | sequencing:RipollSanchez2023 | `data/RipollSanchez2023/NPP_receptor_all_sequencing.csv` |

### 3. Assembly Logic

```
Connectome[source, target] = Release[source, FLP-1] AND Receptor[target, DMSR-7]
where (FLP-1, DMSR-7) is validated as a pair in pairing_info
```

## Intermediate Data Validation

### Pairing Info
- FLP-1 / DMSR-7 pair found in RipollSanchez2023 pairing info
- EC50: 1.42e-11

### Release Vector (FLP-1)
- 2 neurons express FLP-1: **AVKL**, **AVKR**

### Receptor Expression (DMSR-7)
- 102 neurons express DMSR-7
- Includes: motor neurons (DA, VA, DB, VB, AS classes), command interneurons (AVA, AVB, AVD, PVC), ring motor neurons (RMD, SMD), and others

## Results

### Comparison Summary

| Metric | Preassembled | Assembled | Match |
|--------|--------------|-----------|-------|
| Shape | (302, 302) | (302, 302) | Yes |
| Total Connections | 204 | 204 | Yes |
| Source Neurons | 2 | 2 | Yes |
| Target Neurons | 102 | 102 | Yes |
| Exact Matrix Match | - | - | **Yes** |

### Connection Count Verification

```
Expected: 2 sources x 102 targets = 204 connections
Actual: 204 connections
```

### Source Neurons
- Preassembled: AVKL, AVKR
- Assembled: AVKL, AVKR

## Conclusion

**PASSED**: The `assemble_npp_connectome()` function exactly reproduces the preassembled FLP-1/DMSR-7 connectome. All 91,204 matrix entries (302 x 302) match exactly between the assembled and preassembled versions.

## Notes

- The preassembled file uses "longRange" designation, which corresponds to the full neuropeptide signaling range in RipollSanchez2023
- Self-loops are present (AVKL and AVKR express both FLP-1 and DMSR-7)
- The assembled connectome includes connections: AVKL->AVKL, AVKL->AVKR, AVKR->AVKL, AVKR->AVKR
