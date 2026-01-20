"""
Configuration module for molecular connectome assembly.

Provides TOML-based configuration loading, validation, and serialization.
"""

import hashlib
import json
import sys
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Literal

# tomllib is stdlib in Python 3.11+, use tomli for earlier versions
if sys.version_info >= (3, 11):
    import tomllib
else:
    try:
        import tomli as tomllib
    except ImportError:
        raise ImportError(
            "tomli package required for Python < 3.11. "
            "Install with: pip install tomli"
        )


@dataclass
class ReleaseConfig:
    """Configuration for release (neurotransmitter source) data."""
    markers: list[str] = field(default_factory=lambda: ["release"])
    sources: list[str] = field(default_factory=list)
    gate: Literal["or", "and"] = "or"


@dataclass
class ReceptorConfig:
    """Configuration for receptor (target) data."""
    sources: list[str] = field(default_factory=list)
    gate: Literal["or", "and"] = "or"
    type: Literal["all", "ionotropic", "metabotropic"] = "all"


@dataclass
class ConstraintConfig:
    """Configuration for structural constraints."""
    enabled: bool = False
    structural_dataset: str = "Cook2019"
    mode: Literal["binary", "weighted"] = "binary"


@dataclass
class OutputConfig:
    """Configuration for output files."""
    directory: str = "connectomes/candy_assembly"
    basename: str = ""
    formats: list[str] = field(default_factory=lambda: ["csv"])
    save_per_pair: bool = True
    save_count: bool = True
    save_binary: bool = True


@dataclass
class AssemblyConfig:
    """Full configuration for a molecular connectome assembly."""
    version: str = "1.0"
    name: str = ""
    description: str = ""
    author: str = ""
    date: str = ""
    molecule_type: Literal["neurotransmitter", "neuropeptide"] = "neurotransmitter"
    molecule: str = ""
    release: ReleaseConfig = field(default_factory=ReleaseConfig)
    receptor: ReceptorConfig = field(default_factory=ReceptorConfig)
    pairing_source: str | None = None
    constraint: ConstraintConfig = field(default_factory=ConstraintConfig)
    output: OutputConfig = field(default_factory=OutputConfig)


def _parse_release_config(data: dict) -> ReleaseConfig:
    """Parse release section from TOML data."""
    return ReleaseConfig(
        markers=data.get("markers", ["release"]),
        sources=data.get("sources", []),
        gate=data.get("gate", "or"),
    )


def _parse_receptor_config(data: dict) -> ReceptorConfig:
    """Parse receptor section from TOML data."""
    return ReceptorConfig(
        sources=data.get("sources", []),
        gate=data.get("gate", "or"),
        type=data.get("type", "all"),
    )


def _parse_constraint_config(data: dict) -> ConstraintConfig:
    """Parse constraint section from TOML data."""
    return ConstraintConfig(
        enabled=data.get("enabled", False),
        structural_dataset=data.get("structural_dataset", "Cook2019"),
        mode=data.get("mode", "binary"),
    )


def _parse_output_config(data: dict, molecule: str) -> OutputConfig:
    """Parse output section from TOML data."""
    return OutputConfig(
        directory=data.get("directory", f"connectomes/candy_assembly/{molecule}"),
        basename=data.get("basename", molecule),
        formats=data.get("formats", ["csv"]),
        save_per_pair=data.get("save_per_pair", True),
        save_count=data.get("save_count", True),
        save_binary=data.get("save_binary", True),
    )


def load_assembly_config(path: str | Path) -> AssemblyConfig:
    """
    Load an assembly configuration from a TOML file.

    Args:
        path: Path to the TOML configuration file.

    Returns:
        AssemblyConfig dataclass with all configuration fields.

    Raises:
        FileNotFoundError: If the config file doesn't exist.
        ValueError: If required fields are missing.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Config file not found: {path}")

    with open(path, "rb") as f:
        data = tomllib.load(f)

    metadata = data.get("metadata", {})
    assembly = data.get("assembly", {})
    constraint = data.get("constraint", {})
    output = data.get("output", {})

    molecule = assembly.get("molecule", "")

    config = AssemblyConfig(
        version=data.get("version", "1.0"),
        name=metadata.get("name", ""),
        description=metadata.get("description", ""),
        author=metadata.get("author", ""),
        date=metadata.get("date", ""),
        molecule_type=assembly.get("molecule_type", "neurotransmitter"),
        molecule=molecule,
        release=_parse_release_config(assembly.get("release", {})),
        receptor=_parse_receptor_config(assembly.get("receptor", {})),
        pairing_source=assembly.get("pairing_source"),
        constraint=_parse_constraint_config(constraint),
        output=_parse_output_config(output, molecule),
    )

    return config


def validate_config(config: AssemblyConfig) -> list[str]:
    """
    Validate an assembly configuration.

    Args:
        config: AssemblyConfig to validate.

    Returns:
        List of validation error messages. Empty list if valid.
    """
    errors = []

    # Required fields
    if not config.molecule:
        errors.append("assembly.molecule is required")

    if not config.receptor.sources:
        errors.append("assembly.receptor.sources is required (must specify at least one receptor source)")

    # Molecule type validation
    if config.molecule_type not in ("neurotransmitter", "neuropeptide"):
        errors.append(f"assembly.molecule_type must be 'neurotransmitter' or 'neuropeptide', got '{config.molecule_type}'")

    # NPP-specific validation
    if config.molecule_type == "neuropeptide":
        if not config.pairing_source:
            errors.append("assembly.pairing_source is required for neuropeptide assemblies")
        elif config.pairing_source not in ("Altun2013", "Bentley2016", "RipollSanchez2023"):
            errors.append(f"assembly.pairing_source must be one of 'Altun2013', 'Bentley2016', 'RipollSanchez2023', got '{config.pairing_source}'")

    # Gate validation
    if config.release.gate not in ("or", "and"):
        errors.append(f"assembly.release.gate must be 'or' or 'and', got '{config.release.gate}'")

    if config.receptor.gate not in ("or", "and"):
        errors.append(f"assembly.receptor.gate must be 'or' or 'and', got '{config.receptor.gate}'")

    # Receptor type validation
    if config.receptor.type not in ("all", "ionotropic", "metabotropic"):
        errors.append(f"assembly.receptor.type must be 'all', 'ionotropic', or 'metabotropic', got '{config.receptor.type}'")

    # Constraint validation
    if config.constraint.enabled:
        if config.constraint.mode not in ("binary", "weighted"):
            errors.append(f"constraint.mode must be 'binary' or 'weighted', got '{config.constraint.mode}'")

    # Source format validation
    for source in config.release.sources:
        if ":" not in source:
            errors.append(f"Invalid release source format '{source}', expected 'method:dataset'")

    for source in config.receptor.sources:
        if ":" not in source:
            errors.append(f"Invalid receptor source format '{source}', expected 'method:dataset'")

    return errors


def config_to_dict(config: AssemblyConfig) -> dict:
    """Convert AssemblyConfig to a nested dictionary."""
    return asdict(config)


def compute_config_hash(config: AssemblyConfig) -> str:
    """
    Compute a SHA256 hash of the configuration for reproducibility tracking.

    Args:
        config: AssemblyConfig to hash.

    Returns:
        SHA256 hash string prefixed with 'sha256:'.
    """
    config_dict = config_to_dict(config)
    config_json = json.dumps(config_dict, sort_keys=True)
    hash_value = hashlib.sha256(config_json.encode()).hexdigest()
    return f"sha256:{hash_value}"


def save_config(config: AssemblyConfig, path: str | Path) -> None:
    """
    Save an AssemblyConfig to a TOML file.

    Args:
        config: AssemblyConfig to save.
        path: Output path for the TOML file.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    lines = [
        f'version = "{config.version}"',
        "",
        "[metadata]",
        f'name = "{config.name}"',
        f'description = "{config.description}"',
        f'author = "{config.author}"',
        f'date = "{config.date}"',
        "",
        "[assembly]",
        f'molecule_type = "{config.molecule_type}"',
        f'molecule = "{config.molecule}"',
    ]

    if config.pairing_source:
        lines.append(f'pairing_source = "{config.pairing_source}"')

    lines.extend([
        "",
        "[assembly.release]",
        f'markers = {json.dumps(config.release.markers)}',
        f'sources = {json.dumps(config.release.sources)}',
        f'gate = "{config.release.gate}"',
        "",
        "[assembly.receptor]",
        f'sources = {json.dumps(config.receptor.sources)}',
        f'gate = "{config.receptor.gate}"',
        f'type = "{config.receptor.type}"',
        "",
        "[constraint]",
        f'enabled = {"true" if config.constraint.enabled else "false"}',
        f'structural_dataset = "{config.constraint.structural_dataset}"',
        f'mode = "{config.constraint.mode}"',
        "",
        "[output]",
        f'directory = "{config.output.directory}"',
        f'basename = "{config.output.basename}"',
        f'formats = {json.dumps(config.output.formats)}',
        f'save_per_pair = {"true" if config.output.save_per_pair else "false"}',
        f'save_count = {"true" if config.output.save_count else "false"}',
        f'save_binary = {"true" if config.output.save_binary else "false"}',
    ])

    with open(path, "w") as f:
        f.write("\n".join(lines))


def create_metadata_sidecar(
    config: AssemblyConfig,
    output_stats: dict,
    pub_utils_version: str = "0.1.0",
) -> dict:
    """
    Create a metadata sidecar dictionary for an assembly output.

    Args:
        config: AssemblyConfig used for the assembly.
        output_stats: Statistics from the assembly (e.g., total_connections, receptors_included).
        pub_utils_version: Version of pub_utils used.

    Returns:
        Dictionary containing all metadata for the assembly.
    """
    return {
        "assembly_name": config.name,
        "timestamp": datetime.now().isoformat(),
        "config_hash": compute_config_hash(config),
        "pub_utils_version": pub_utils_version,
        "config": config_to_dict(config),
        "output_stats": output_stats,
    }
