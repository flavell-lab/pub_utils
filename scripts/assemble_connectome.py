#!/usr/bin/env python3
"""
CLI runner for molecular connectome assembly.

Usage:
    python scripts/assemble_connectome.py configs/dopamine_reporter.toml

This script:
1. Loads and validates a TOML configuration file
2. Calls the appropriate assembly function (NT or NPP)
3. Applies structural constraints if enabled
4. Saves outputs with a metadata sidecar
"""

import argparse
import json
import sys
import warnings
from pathlib import Path

# Add src to path for development use
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

import pub_utils as pu
from pub_utils.config import (
    load_assembly_config,
    validate_config,
    create_metadata_sidecar,
)


def run_assembly(config_path: str, verbose: bool = True) -> dict:
    """
    Run a molecular connectome assembly from a TOML configuration file.

    Args:
        config_path: Path to the TOML configuration file.
        verbose: Whether to print progress messages.

    Returns:
        dict with keys:
            'assembly': the assembly result dict
            'config': the AssemblyConfig used
            'output_files': list of output file paths created
    """
    # Load and validate config
    if verbose:
        print(f"Loading config: {config_path}")

    config = load_assembly_config(config_path)
    errors = validate_config(config)

    if errors:
        print("Configuration validation errors:")
        for error in errors:
            print(f"  - {error}")
        sys.exit(1)

    if verbose:
        print(f"Assembling {config.molecule_type}: {config.molecule}")
        print(f"  Release sources: {config.release.sources or 'all available'}")
        print(f"  Receptor sources: {config.receptor.sources}")

    # Run assembly
    if config.molecule_type == "neurotransmitter":
        assembly = pu.assemble_nt_connectome(
            neurotransmitter=config.molecule,
            release_markers=config.release.markers,
            release_sources=config.release.sources or None,
            release_gate=config.release.gate,
            receptor_sources=config.receptor.sources,
            receptor_gate=config.receptor.gate,
            receptor_type=config.receptor.type,
        )
    else:  # neuropeptide
        assembly = pu.assemble_npp_connectome(
            neuropeptide=config.molecule,
            release_sources=config.release.sources or None,
            release_gate=config.release.gate,
            receptor_sources=config.receptor.sources,
            pairing_source=config.pairing_source,
            receptor_gate=config.receptor.gate,
        )

    # Apply structural constraints if enabled
    if config.constraint.enabled:
        if verbose:
            print(f"  Applying structural constraint: {config.constraint.structural_dataset} ({config.constraint.mode})")

        assembly = pu.constrain_assembly(
            assembly,
            structural_dataset=config.constraint.structural_dataset,
            mode=config.constraint.mode,
        )

    # Compute output statistics
    binary_matrix = assembly["binary"]
    count_matrix = assembly["count"]

    total_connections = int((binary_matrix == 1).sum().sum())
    receptors_included = list(assembly["per_pair"].keys())

    output_stats = {
        "total_connections": total_connections,
        "receptors_included": receptors_included,
        "matrix_shape": list(binary_matrix.shape),
    }

    if verbose:
        print(f"  Total connections: {total_connections}")
        print(f"  Receptors included: {len(receptors_included)}")

    # Save outputs
    output_dir = Path(config.output.directory)
    output_dir.mkdir(parents=True, exist_ok=True)

    basename = config.output.basename or config.molecule
    output_files = []

    if config.output.save_binary:
        binary_path = output_dir / f"{basename}_binary.csv"
        assembly["binary"].to_csv(binary_path)
        output_files.append(str(binary_path))
        if verbose:
            print(f"  Saved: {binary_path}")

    if config.output.save_count:
        count_path = output_dir / f"{basename}_count.csv"
        assembly["count"].to_csv(count_path)
        output_files.append(str(count_path))
        if verbose:
            print(f"  Saved: {count_path}")

    if config.output.save_per_pair:
        per_pair_dir = output_dir / f"{basename}_per_pair"
        per_pair_dir.mkdir(parents=True, exist_ok=True)
        for receptor_name, matrix in assembly["per_pair"].items():
            receptor_path = per_pair_dir / f"{receptor_name}.csv"
            matrix.to_csv(receptor_path)
            output_files.append(str(receptor_path))
        if verbose:
            print(f"  Saved per-pair matrices to: {per_pair_dir}/")

    # Save metadata sidecar
    metadata = create_metadata_sidecar(config, output_stats, pu.__version__)
    metadata_path = output_dir / f"{basename}_metadata.json"
    with open(metadata_path, "w") as f:
        json.dump(metadata, f, indent=2)
    output_files.append(str(metadata_path))
    if verbose:
        print(f"  Saved: {metadata_path}")

    return {
        "assembly": assembly,
        "config": config,
        "output_files": output_files,
    }


def main():
    parser = argparse.ArgumentParser(
        description="Assemble a molecular connectome from a TOML configuration file."
    )
    parser.add_argument(
        "config",
        type=str,
        help="Path to the TOML configuration file",
    )
    parser.add_argument(
        "-q", "--quiet",
        action="store_true",
        help="Suppress progress messages",
    )
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Only validate the config, don't run the assembly",
    )

    args = parser.parse_args()

    if args.validate_only:
        config = load_assembly_config(args.config)
        errors = validate_config(config)
        if errors:
            print("Validation errors:")
            for error in errors:
                print(f"  - {error}")
            sys.exit(1)
        else:
            print("Configuration is valid.")
            sys.exit(0)

    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        result = run_assembly(args.config, verbose=not args.quiet)

        # Print any warnings that occurred
        if w and not args.quiet:
            print("\nWarnings during assembly:")
            for warning in w:
                print(f"  - {warning.message}")

    if not args.quiet:
        print("\nAssembly complete.")


if __name__ == "__main__":
    main()
