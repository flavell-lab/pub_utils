# Molecular Connectome Assembly with Config-Driven Workflow

## Overview

This document describes a standardized, reproducible approach for assembling molecular connectomes using the existing `pub_utils` infrastructure. The recommended approach is **TOML configuration files with a runner script**, which provides declarative documentation, version control, and deterministic execution.

## Why Config-Driven Assembly?

| Benefit | Description |
|---------|-------------|
| Reproducibility | Configs are version-controlled and self-documenting |
| Determinism | Same config always produces identical output |
| Claude Code Integration | Claude can generate, validate, and execute configs |
| Low Complexity | Uses existing Python functions without new infrastructure |

---

## Implementation Phases

### Phase 1: Config Module (`src/pub_utils/config.py`)

Create a configuration schema and validation module using TOML format.

**Key functions to implement:**

```python
from dataclasses import dataclass
from typing import Literal

@dataclass
class ReleaseConfig:
    markers: list[str]           # Functional markers or gene names
    sources: list[str]           # method:dataset format (e.g., "reporter:Wang2024")
    gate: Literal["or", "and"]   # Logical gate for combining sources

@dataclass
class ReceptorConfig:
    sources: list[str]           # method:dataset format
    gate: Literal["or", "and"]
    type: Literal["all", "ionotropic", "metabotropic"]

@dataclass
class ConstraintConfig:
    enabled: bool
    structural_dataset: str      # e.g., "Cook2019"
    mode: Literal["binary", "weighted"]

@dataclass
class OutputConfig:
    directory: str
    basename: str
    formats: list[str]           # e.g., ["csv"]
    save_per_pair: bool
    save_count: bool
    save_binary: bool

@dataclass
class AssemblyConfig:
    version: str
    name: str
    description: str
    author: str
    date: str
    molecule_type: Literal["neurotransmitter", "neuropeptide"]
    molecule: str
    release: ReleaseConfig
    receptor: ReceptorConfig
    pairing_source: str | None   # For NPP: 'Altun2013', 'Bentley2016', 'RipollSanchez2023'
    constraint: ConstraintConfig
    output: OutputConfig

def load_assembly_config(path: str) -> AssemblyConfig
def validate_config(config: AssemblyConfig) -> list[str]  # Returns validation errors
def save_config(config: AssemblyConfig, path: str) -> None
```

**TOML Config Schema:**

```toml
version = "1.0"

[metadata]
name = "dopamine_reporter_constrained"
description = "Dopamine connectome using Wang2024 reporter data"
author = "researcher_name"
date = "2026-01-19"

[assembly]
molecule_type = "neurotransmitter"  # or "neuropeptide"
molecule = "dopamine"

[assembly.release]
markers = ["release"]               # Functional markers or gene names
sources = ["reporter:Wang2024"]     # method:dataset format
gate = "or"                         # 'or' or 'and'

[assembly.receptor]
sources = ["reporter:Muralidhara2025"]
gate = "or"
type = "all"                        # 'all', 'ionotropic', 'metabotropic'

# pairing_source: For NPP only - 'Altun2013', 'Bentley2016', 'RipollSanchez2023'
# pairing_source = "RipollSanchez2023"

[constraint]
enabled = true
structural_dataset = "Cook2019"
mode = "binary"                     # 'binary' or 'weighted'

[output]
directory = "connectomes/candy_assembly/dopamine"
basename = "dopamine_reporter_Cook2019"
formats = ["csv"]
save_per_pair = true
save_count = true
save_binary = true
```

---

### Phase 2: Runner Script (`scripts/assemble_connectome.py`)

Create a CLI runner that:
1. Loads and validates TOML config
2. Calls existing `assemble_nt_connectome()` or `assemble_npp_connectome()`
3. Applies structural constraints if enabled
4. Saves outputs with metadata sidecar

**Usage:**
```bash
python scripts/assemble_connectome.py configs/dopamine_reporter.toml
```

**Metadata Sidecar (`{basename}_metadata.json`):**
```json
{
  "assembly_name": "dopamine_reporter_Cook2019",
  "timestamp": "2026-01-19T14:30:00Z",
  "config_hash": "sha256:abc123...",
  "pub_utils_version": "0.1.0",
  "config": { /* embedded config */ },
  "output_stats": {
    "total_connections": 156,
    "receptors_included": ["dop-1", "dop-2", "dop-3", "dop-4", "dop-5", "dop-6"]
  }
}
```

---

### Phase 3: Example Configs and Documentation

Create example configs in `configs/examples/`:
- `nt_dopamine_example.toml`
- `nt_acetylcholine_ionotropic_example.toml`
- `npp_flp1_example.toml`

---

### Phase 4: MCP Server (Future Implementation)

An MCP server will enable direct Claude Code tool integration for assembly operations.

**File:** `src/pub_utils/mcp_server.py`

**Tools to expose:**

| Tool | Description |
|------|-------------|
| `list_available_sources` | List all data sources from assets.json |
| `validate_assembly_config` | Validate a TOML config |
| `assemble_connectome` | Run assembly from config |
| `compare_assemblies` | Compare two connectome outputs |

**MCP Server Implementation Outline:**

```python
from mcp.server import Server
from mcp.server.stdio import stdio_server
from mcp.types import Tool, TextContent

server = Server("pub_utils")

@server.list_tools()
async def list_tools() -> list[Tool]:
    return [
        Tool(
            name="list_available_sources",
            description="List all available data sources from assets.json",
            inputSchema={
                "type": "object",
                "properties": {
                    "molecule_type": {
                        "type": "string",
                        "enum": ["neurotransmitter", "neuropeptide"],
                        "description": "Type of molecule to list sources for"
                    }
                }
            }
        ),
        Tool(
            name="validate_assembly_config",
            description="Validate a TOML assembly configuration",
            inputSchema={
                "type": "object",
                "properties": {
                    "config_path": {
                        "type": "string",
                        "description": "Path to TOML config file"
                    }
                },
                "required": ["config_path"]
            }
        ),
        Tool(
            name="assemble_connectome",
            description="Assemble a molecular connectome from a validated config",
            inputSchema={
                "type": "object",
                "properties": {
                    "config_path": {
                        "type": "string",
                        "description": "Path to TOML config file"
                    }
                },
                "required": ["config_path"]
            }
        ),
        Tool(
            name="compare_assemblies",
            description="Compare two connectome outputs",
            inputSchema={
                "type": "object",
                "properties": {
                    "connectome_a": {"type": "string"},
                    "connectome_b": {"type": "string"}
                },
                "required": ["connectome_a", "connectome_b"]
            }
        )
    ]

@server.call_tool()
async def call_tool(name: str, arguments: dict) -> list[TextContent]:
    # Implementation for each tool
    pass

async def main():
    async with stdio_server() as (read_stream, write_stream):
        await server.run(read_stream, write_stream)
```

**MCP Configuration (for Claude Desktop / Claude Code):**

```json
{
  "mcpServers": {
    "pub_utils": {
      "command": "python",
      "args": ["-m", "pub_utils.mcp_server"],
      "cwd": "/path/to/pub_utils"
    }
  }
}
```

---

## Critical Files to Create/Modify

| File | Action | Phase |
|------|--------|-------|
| `src/pub_utils/config.py` | Create new - config loading/validation | 1 |
| `src/pub_utils/__init__.py` | Add exports for config functions | 1 |
| `scripts/assemble_connectome.py` | Create new - CLI runner | 2 |
| `configs/examples/*.toml` | Create new - example configs | 3 |
| `src/pub_utils/mcp_server.py` | Create new - MCP server | 4 (future) |

## Existing Functions to Use (No Modification)

- `pub_utils.assemble_nt_connectome()` - Core NT assembly
- `pub_utils.assemble_npp_connectome()` - Core NPP assembly
- `pub_utils.constrain_assembly()` - Structural constraint application
- `pub_utils.constants.AllHermNeurons` - Canonical neuron ordering

---

## Verification Plan

### Unit Tests (`tests/test_config.py`)
- Config schema validation
- Invalid config rejection
- Source string parsing (method:dataset format)

### Integration Tests (`tests/test_assembly_runner.py`)
- Round-trip: config → assembly → validate
- Compare to preassembled reference connectomes (e.g., `dopamine_Bentley2016.csv`)
- Determinism: same config produces identical output

### Validation Checklist for Each Assembly
1. Connection count within expected range
2. Matrix shape is 302×302 (hermaphrodite neurons)
3. Comparison to similar assemblies shows reasonable overlap
4. Metadata sidecar is complete and accurate

---

## Example Workflow with Claude Code

```bash
# 1. Claude generates config based on user requirements
# User: "Assemble a dopamine connectome using Wang2024 reporter data"

# 2. Config file created at configs/dopamine_wang2024.toml

# 3. Run assembly
python scripts/assemble_connectome.py configs/dopamine_wang2024.toml

# 4. Output files:
#    - connectomes/candy_assembly/dopamine/dopamine_wang2024_binary.csv
#    - connectomes/candy_assembly/dopamine/dopamine_wang2024_count.csv
#    - connectomes/candy_assembly/dopamine/dopamine_wang2024_metadata.json

# 5. Log to claudecode/2026-01-19.md
```

---

## Summary

| Phase | Deliverable | Status |
|-------|-------------|--------|
| 1 | Config module (`src/pub_utils/config.py`) | To implement |
| 2 | Runner script (`scripts/assemble_connectome.py`) | To implement |
| 3 | Example TOML configs | To implement |
| 4 | MCP Server | Future |

The config-driven approach ensures reproducibility through declarative, version-controlled assembly specifications while leveraging the existing `pub_utils` infrastructure.
