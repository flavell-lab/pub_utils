#!/usr/bin/env python3
"""
MCP Server for C. elegans Neuron Similarity

Exposes neuron similarity computation and data access tools via Model Context Protocol.

Usage:
    python mcp_server/neuron_similarity.py

Configuration (add to ~/.claude.json or .claude/settings.json):
    {
      "mcpServers": {
        "neuron-similarity": {
          "command": "python",
          "args": ["/path/to/pub_utils/mcp_server/neuron_similarity.py"]
        }
      }
    }
"""

import json
import sys
from pathlib import Path

# Add src directory to path for imports
_REPO_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(_REPO_ROOT / "src"))

from mcp.server import Server
from mcp.server.stdio import stdio_server
from mcp.types import Tool, TextContent, Resource

import pandas as pd

# Import pub_utils modules
import pub_utils as pu
from pub_utils.similarity import (
    get_sensory_profile,
    get_process_profile,
    get_nt_release_profile,
    get_nt_receptor_profile,
    get_npp_release_profile,
    get_npp_receptor_profile,
    compute_sensory_similarity,
    compute_process_similarity,
    compute_release_similarity,
    compute_receptor_similarity,
    find_similar_neurons,
    compare_similarity_dimensions,
)
from pub_utils.constants import AllHermNeurons

# Initialize server
server = Server("neuron-similarity")

# Cache for computed similarity matrices
_SIMILARITY_CACHE = {}


# =============================================================================
# Tools
# =============================================================================

@server.tool()
async def get_neuron_features(neuron_id: str) -> str:
    """
    Get all features for a specific neuron.

    Args:
        neuron_id: The neuron identifier (e.g., 'ASEL', 'AVAR', 'DA01')

    Returns:
        JSON with cell type, sensory types, segment, and anatomical processes
    """
    if neuron_id not in AllHermNeurons:
        return json.dumps({"error": f"Neuron '{neuron_id}' not found. Use one of the 302 hermaphrodite neurons."})

    sensory = get_sensory_profile()
    process = get_process_profile()

    sensory_types = sensory.loc[neuron_id]
    active_sensory = sensory_types[sensory_types == 1].index.tolist()

    process_locs = process.loc[neuron_id]
    active_processes = process_locs[process_locs == 1].index.tolist()

    return json.dumps({
        "neuron_id": neuron_id,
        "sensory_types": active_sensory,
        "anatomical_processes": active_processes,
        "n_sensory_types": len(active_sensory),
        "n_processes": len(active_processes)
    }, indent=2)


@server.tool()
async def get_release_profile(
    neuron_id: str,
    molecule_type: str = "all"
) -> str:
    """
    Get neurotransmitter/neuropeptide release profile for a neuron.

    Args:
        neuron_id: The neuron identifier
        molecule_type: 'nt' (neurotransmitters), 'npp' (neuropeptides), or 'all'

    Returns:
        JSON with released molecules
    """
    if neuron_id not in AllHermNeurons:
        return json.dumps({"error": f"Neuron '{neuron_id}' not found"})

    result = {"neuron_id": neuron_id}

    if molecule_type in ["nt", "all"]:
        nt_profile = get_nt_release_profile()
        nt_row = nt_profile.loc[neuron_id]
        released_nt = nt_row[nt_row == 1].index.tolist()
        result["neurotransmitters"] = released_nt

    if molecule_type in ["npp", "all"]:
        npp_profile = get_npp_release_profile()
        npp_row = npp_profile.loc[neuron_id]
        released_npp = npp_row[npp_row == 1].index.tolist()
        result["neuropeptides"] = released_npp

    return json.dumps(result, indent=2)


@server.tool()
async def get_receptor_profile(
    neuron_id: str,
    molecule_type: str = "all"
) -> str:
    """
    Get receptor expression profile for a neuron.

    Args:
        neuron_id: The neuron identifier
        molecule_type: 'nt' (neurotransmitter receptors), 'npp' (neuropeptide receptors), or 'all'

    Returns:
        JSON with expressed receptors
    """
    if neuron_id not in AllHermNeurons:
        return json.dumps({"error": f"Neuron '{neuron_id}' not found"})

    result = {"neuron_id": neuron_id}

    if molecule_type in ["nt", "all"]:
        nt_profile = get_nt_receptor_profile()
        if not nt_profile.empty:
            nt_row = nt_profile.loc[neuron_id]
            expressed_nt = nt_row[nt_row == 1].index.tolist()
            result["nt_receptors"] = expressed_nt

    if molecule_type in ["npp", "all"]:
        npp_profile = get_npp_receptor_profile()
        if not npp_profile.empty:
            npp_row = npp_profile.loc[neuron_id]
            expressed_npp = npp_row[npp_row == 1].index.tolist()
            result["npp_receptors"] = expressed_npp

    return json.dumps(result, indent=2)


@server.tool()
async def compute_similarity(
    neuron_a: str,
    neuron_b: str,
    dimension: str = "all"
) -> str:
    """
    Compute similarity between two neurons.

    Args:
        neuron_a: First neuron identifier
        neuron_b: Second neuron identifier
        dimension: 'sensory', 'process', 'release', 'receptor', or 'all'

    Returns:
        JSON with similarity scores
    """
    for n in [neuron_a, neuron_b]:
        if n not in AllHermNeurons:
            return json.dumps({"error": f"Neuron '{n}' not found"})

    results = {}

    dimensions_to_compute = [dimension] if dimension != "all" else ["sensory", "process", "release", "receptor"]

    for dim in dimensions_to_compute:
        # Check cache
        if dim not in _SIMILARITY_CACHE:
            if dim == "sensory":
                _SIMILARITY_CACHE[dim] = compute_sensory_similarity()
            elif dim == "process":
                _SIMILARITY_CACHE[dim] = compute_process_similarity()
            elif dim == "release":
                _SIMILARITY_CACHE[dim] = compute_release_similarity()
            elif dim == "receptor":
                _SIMILARITY_CACHE[dim] = compute_receptor_similarity()

        sim_matrix = _SIMILARITY_CACHE[dim]["similarity"]
        results[dim] = float(sim_matrix.loc[neuron_a, neuron_b])

    return json.dumps({
        "neuron_a": neuron_a,
        "neuron_b": neuron_b,
        "similarity": results
    }, indent=2)


@server.tool()
async def find_similar(
    neuron_id: str,
    dimension: str = "sensory",
    top_n: int = 10
) -> str:
    """
    Find neurons most similar to a query neuron.

    Args:
        neuron_id: The query neuron identifier
        dimension: 'sensory', 'process', 'release', or 'receptor'
        top_n: Number of similar neurons to return (max 50)

    Returns:
        JSON with list of similar neurons and their similarity scores
    """
    if neuron_id not in AllHermNeurons:
        return json.dumps({"error": f"Neuron '{neuron_id}' not found"})

    top_n = min(top_n, 50)

    # Check cache
    if dimension not in _SIMILARITY_CACHE:
        if dimension == "sensory":
            _SIMILARITY_CACHE[dimension] = compute_sensory_similarity()
        elif dimension == "process":
            _SIMILARITY_CACHE[dimension] = compute_process_similarity()
        elif dimension == "release":
            _SIMILARITY_CACHE[dimension] = compute_release_similarity()
        elif dimension == "receptor":
            _SIMILARITY_CACHE[dimension] = compute_receptor_similarity()
        else:
            return json.dumps({"error": f"Unknown dimension '{dimension}'"})

    similar_df = find_similar_neurons(
        _SIMILARITY_CACHE[dimension]["similarity"],
        neuron_id,
        top_n=top_n
    )

    return json.dumps({
        "query_neuron": neuron_id,
        "dimension": dimension,
        "similar_neurons": [
            {"neuron_id": row["neuron_id"], "similarity": round(row["similarity"], 4)}
            for _, row in similar_df.iterrows()
        ]
    }, indent=2)


@server.tool()
async def list_available_data() -> str:
    """
    List available data sources and their coverage.

    Returns:
        JSON with available data sources, neuron counts, and feature counts
    """
    sensory = get_sensory_profile()
    process = get_process_profile()

    return json.dumps({
        "neurons": {
            "total": len(AllHermNeurons),
            "list_sample": AllHermNeurons[:10] + ["..."]
        },
        "features": {
            "sensory_types": {
                "count": sensory.shape[1],
                "names": sensory.columns.tolist()
            },
            "anatomical_processes": {
                "count": process.shape[1],
                "names": process.columns.tolist()
            }
        },
        "similarity_dimensions": [
            "sensory - similarity based on sensory modality responses",
            "process - similarity based on anatomical process locations",
            "release - similarity based on NT/NPP release profiles",
            "receptor - similarity based on receptor expression"
        ]
    }, indent=2)


@server.tool()
async def get_neurons_by_sensory_type(sensory_type: str) -> str:
    """
    Get all neurons that respond to a specific sensory modality.

    Args:
        sensory_type: One of: chemical, odor, pheromone, osmolarity, oxygen,
                     carbon_dioxide, photon, thermal, electrical, mechanical,
                     stretch, nociceptive

    Returns:
        JSON with list of neurons
    """
    sensory = get_sensory_profile()

    if sensory_type not in sensory.columns:
        return json.dumps({
            "error": f"Unknown sensory type '{sensory_type}'",
            "available_types": sensory.columns.tolist()
        })

    neurons = sensory[sensory[sensory_type] == 1].index.tolist()

    return json.dumps({
        "sensory_type": sensory_type,
        "neurons": neurons,
        "count": len(neurons)
    }, indent=2)


@server.tool()
async def get_neurons_by_process(process_name: str) -> str:
    """
    Get all neurons that have processes in a specific anatomical location.

    Args:
        process_name: Anatomical location (e.g., 'nerve_ring', 'dorsal_cord')

    Returns:
        JSON with list of neurons
    """
    process = get_process_profile()

    if process_name not in process.columns:
        return json.dumps({
            "error": f"Unknown process location '{process_name}'",
            "available_locations": process.columns.tolist()
        })

    neurons = process[process[process_name] == 1].index.tolist()

    return json.dumps({
        "process_location": process_name,
        "neurons": neurons,
        "count": len(neurons)
    }, indent=2)


# =============================================================================
# Resources
# =============================================================================

@server.resource("similarity://sensory")
async def get_sensory_similarity_matrix() -> str:
    """Get the full sensory similarity matrix as CSV."""
    if "sensory" not in _SIMILARITY_CACHE:
        _SIMILARITY_CACHE["sensory"] = compute_sensory_similarity()
    return _SIMILARITY_CACHE["sensory"]["similarity"].to_csv()


@server.resource("similarity://process")
async def get_process_similarity_matrix() -> str:
    """Get the full anatomical process similarity matrix as CSV."""
    if "process" not in _SIMILARITY_CACHE:
        _SIMILARITY_CACHE["process"] = compute_process_similarity()
    return _SIMILARITY_CACHE["process"]["similarity"].to_csv()


@server.resource("neurons://list")
async def get_neuron_list() -> str:
    """Get the list of all 302 hermaphrodite neurons."""
    return json.dumps(AllHermNeurons)


@server.resource("features://sensory")
async def get_sensory_features() -> str:
    """Get available sensory type features."""
    sensory = get_sensory_profile()
    return json.dumps({
        "features": sensory.columns.tolist(),
        "count": len(sensory.columns)
    })


@server.resource("features://process")
async def get_process_features() -> str:
    """Get available anatomical process features."""
    process = get_process_profile()
    return json.dumps({
        "features": process.columns.tolist(),
        "count": len(process.columns)
    })


# =============================================================================
# Main
# =============================================================================

async def main():
    """Run the MCP server."""
    async with stdio_server() as (read_stream, write_stream):
        await server.run(read_stream, write_stream, server.create_initialization_options())


if __name__ == "__main__":
    import asyncio
    asyncio.run(main())
