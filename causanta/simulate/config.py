"""Configuration loading and validation for CAUSANTA simulation."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class DomainConfig:
    width_um: int = 1000
    height_um: int = 1000
    env_grid_um: int = 10


@dataclass(frozen=True)
class TimeConfig:
    total_hours: int = 168
    dt_diffusion_hr: float = 0.01
    output_interval_hr: int = 1
    burnin_hours: int = 100


@dataclass(frozen=True)
class SubstrateConfig:
    diffusion_coeff_um2_hr: float = 0.0
    decay_rate_per_hr: float = 0.0
    boundary_value: float = 0.0
    initial_value: float = 0.0
    clamp_min: float = 0.0
    clamp_max: float = 100.0


@dataclass(frozen=True)
class EnvironmentConfig:
    O2_blood_mmHg: float = 60.0
    q_O2_transfer_per_hr: float = 5.0
    q_glucose_transfer_per_hr: float = 3.0
    hypoxia_threshold_mmHg: float = 10.0
    substrates: dict[str, SubstrateConfig] = field(default_factory=dict)
    # Hypoxia effect on EGFR expression (HIF-2alpha mediated translation).
    # Set to 0.0 to remove the Hypoxia→EGFR edge for robustness analysis.
    # Reference: Franovic et al. 2007 PNAS showed ~2-3x increase.
    egfr_hypoxia_upregulation: float = 1.5  # (1+value)x increase under hypoxia


@dataclass(frozen=True)
class AngiogenesisConfig:
    angiogenesis_threshold_nM: float = 5.0
    vessel_maturation_hr: float = 48.0
    sprouting_rate_per_hr: float = 0.01
    max_sprout_length: int = 200


@dataclass(frozen=True)
class ImmuneRecruitmentConfig:
    """Configuration for immune cell recruitment and killing dynamics.

    Based on experimental literature:
    - Immunological synapse formation: 1-2 hours (Huppa & Davis 2003, Dustin 2008)
    - CTL serial killing capacity: ~10 cells before exhaustion (Boissonnas 2007)
    - Chemotaxis sensitivity: follows CCL2/CCL5/CXCL10 gradients
    - Kill probability: ~70% per successful synapse (varies with activation)
    """

    # Recruitment parameters
    chemokine_threshold: float = 0.5  # Threshold for immune cell recruitment
    recruitment_rate_per_hr: float = 0.005  # Rate of new immune cell entry

    # Killing mechanics
    kill_radius_um: float = 20.0  # Distance for detecting kill targets
    synapse_formation_time_hr: float = 1.5  # Time to form immunological synapse (1-2 hrs)
    kill_probability_per_synapse: float = 0.7  # P(kill) once synapse forms
    min_contact_for_kill_hr: float = 0.5  # Minimum contact before any kill possible

    # Exhaustion mechanics (CTL serial killing limits)
    max_kills_before_exhaustion: int = 10  # CTLs exhaust after ~10 kills
    exhaustion_per_kill: float = 0.1  # Exhaustion increment per successful kill
    exhaustion_recovery_rate_per_hr: float = 0.02  # Slow recovery when not killing
    exhausted_kill_penalty: float = 0.5  # Kill probability multiplier when exhausted

    # Activation dynamics
    activation_radius_um: float = 50.0  # Distance to detect tumor for activation
    activation_rate_per_hr: float = 0.3  # Rate of activation when near tumor
    deactivation_rate_per_hr: float = 0.05  # Rate of deactivation when away from tumor
    min_activation_for_kill: float = 0.2  # Minimum activation level to attempt kill

    # Chemotaxis parameters
    chemokine_diffusion_um2_hr: float = 1800.0  # ~0.5 um²/s typical for chemokines
    chemokine_decay_per_hr: float = 0.5  # Half-life ~1.4 hours
    tumor_chemokine_secretion: float = 1.0  # Relative secretion rate from tumor cells


@dataclass(frozen=True)
class TransitionRule:
    condition: str = ""
    target_type: int = -1
    subtype_flag: str = ""
    rate_per_hr: float = 0.0


@dataclass(frozen=True)
class CellTypeConfig:
    type_id: int = 0
    type_name: str = ""
    shape_path: str = ""
    color: str = "#FFFFFF"
    diameter_mean_um: float = 10.0
    diameter_std_um: float = 2.0
    can_divide: bool = False
    division_time_mean_hr: float = 0.0
    division_time_std_hr: float = 0.0
    O2_prolif_threshold_mmHg: float = 0.0
    glucose_prolif_threshold_mM: float = 0.0
    max_generations: int = -1
    migration_speed_um_hr: float = 0.0
    migration_persistence_hr: float = 1.0
    chemotaxis_O2: float = 0.0
    chemotaxis_VEGF: float = 0.0
    chemotaxis_chemokine: float = 0.0
    haptotaxis_ECM: float = 0.0
    contact_inhibited: bool = False
    O2_necrosis_threshold_mmHg: float = 2.5
    necrosis_delay_hr: float = 6.0
    apoptosis_rate_per_hr: float = 1e-4
    lysis_rate_per_hr: float = 0.01
    O2_consumption_amol_hr: float = 0.0
    glucose_consumption_amol_hr: float = 0.0
    VEGF_secretion_amol_hr: float = 0.0
    lactate_production_amol_hr: float = 0.0
    transition_rules: tuple[TransitionRule, ...] = ()
    ecDNA_init_count: int = 0
    ecDNA_cargo: tuple[str, ...] = ()
    ecDNA_segregation_p: float = 0.5
    ecDNA_effect_on_division: float = 0.0
    ecDNA_effect_on_VEGF: float = 0.0
    ecDNA_effect_on_migration: float = 0.0
    ecDNA_effect_on_survival: float = 0.0


@dataclass(frozen=True)
class InitializationConfig:
    cell_density_per_um2: float = 0.005
    fractions: dict[str, float] = field(default_factory=dict)
    vascular_spacing_um: float = 150.0
    astrocyte_domain_um: float = 75.0


@dataclass(frozen=True)
class TumorSeed:
    x: int = 500
    y: int = 500
    ecDNA_count: int = 20
    ecDNA_cargo: tuple[str, ...] = ("EGFR", "MYC")
    cell_cycle_phase: str = "G1"


@dataclass(frozen=True)
class SimulationConfig:
    rng_seed: int = 42
    output_dir: str | None = None  # Explicit output directory; if None, timestamped dir is created
    domain: DomainConfig = field(default_factory=DomainConfig)
    time: TimeConfig = field(default_factory=TimeConfig)
    environment: EnvironmentConfig = field(default_factory=EnvironmentConfig)
    angiogenesis: AngiogenesisConfig = field(default_factory=AngiogenesisConfig)
    immune_recruitment: ImmuneRecruitmentConfig = field(default_factory=ImmuneRecruitmentConfig)
    cell_types: dict[int, CellTypeConfig] = field(default_factory=dict)
    tumor_seeds: tuple[TumorSeed, ...] = ()
    initialization: InitializationConfig = field(default_factory=InitializationConfig)


# --- Shape paths per cell type (spec Section 3.3) ---
DEFAULT_SHAPE_PATHS = {
    "Neuron": "M 0 -1 C 0.5 -0.5 0.3 0.2 0.8 1 L 0 0.6 L -0.8 1 C -0.3 0.2 -0.5 -0.5 0 -1 Z",
    "Astrocyte": "M 0 -1 L 0.4 -0.3 L 1 -0.2 L 0.4 0.2 L 0.5 1 L 0 0.5 L -0.5 1 L -0.4 0.2 L -1 -0.2 L -0.4 -0.3 Z",
    "Oligodendrocyte": "M 0 -0.8 A 0.8 0.8 0 1 1 0 0.8 A 0.8 0.8 0 1 1 0 -0.8 Z",
    "Microglia": "M 0 -0.6 L 0.3 -0.3 L 0.8 -0.5 L 0.5 0 L 0.8 0.5 L 0.3 0.3 L 0 0.6 L -0.3 0.3 L -0.8 0.5 L -0.5 0 L -0.8 -0.5 L -0.3 -0.3 Z",
    "Endothelial": "M -1 -0.3 Q 0 -0.5 1 -0.3 L 1 0.3 Q 0 0.5 -1 0.3 Z",
    "Pericyte": "M -0.6 -0.4 Q 0 -0.6 0.6 -0.4 L 0.6 0.4 Q 0 0.6 -0.6 0.4 Z",
    "Tumor": "M 0 -0.9 C 0.6 -0.9 0.9 -0.4 0.9 0 C 0.9 0.5 0.5 0.9 0 0.9 C -0.5 0.9 -0.9 0.5 -0.9 0 C -0.9 -0.4 -0.6 -0.9 0 -0.9 Z",
    "RecruitedImmune": "M 0 -0.7 A 0.7 0.7 0 1 1 0 0.7 A 0.7 0.7 0 1 1 0 -0.7 Z",
    "Necrotic": "M -0.5 -0.5 L 0.5 -0.5 L 0.5 0.5 L -0.5 0.5 Z",
}


def _parse_substrate(data: dict[str, Any]) -> SubstrateConfig:
    return SubstrateConfig(**{k: data[k] for k in SubstrateConfig.__dataclass_fields__ if k in data})


def _parse_transition_rules(rules_data: list[dict]) -> tuple[TransitionRule, ...]:
    result = []
    for r in rules_data:
        result.append(TransitionRule(
            condition=r.get("condition", ""),
            target_type=r.get("target_type", -1),
            subtype_flag=r.get("subtype_flag", ""),
            rate_per_hr=r.get("rate_per_hr", 0.0),
        ))
    return tuple(result)


def _parse_cell_type(type_id: int, data: dict[str, Any]) -> CellTypeConfig:
    kwargs: dict[str, Any] = {"type_id": type_id}
    for k, f in CellTypeConfig.__dataclass_fields__.items():
        if k == "type_id":
            continue
        if k == "transition_rules":
            kwargs[k] = _parse_transition_rules(data.get(k, []))
        elif k == "ecDNA_cargo":
            kwargs[k] = tuple(data.get(k, []))
        elif k in data:
            kwargs[k] = data[k]
    return CellTypeConfig(**kwargs)


def _parse_tumor_seeds(seeds_data: list[dict]) -> tuple[TumorSeed, ...]:
    result = []
    for s in seeds_data:
        result.append(TumorSeed(
            x=s.get("x", 500),
            y=s.get("y", 500),
            ecDNA_count=s.get("ecDNA_count", 20),
            ecDNA_cargo=tuple(s.get("ecDNA_cargo", ["EGFR", "MYC"])),
            cell_cycle_phase=s.get("cell_cycle_phase", "G1"),
        ))
    return tuple(result)


def load_config(json_path: str | Path) -> SimulationConfig:
    """Load simulation parameters from JSON file, validate, return frozen config."""
    path = Path(json_path)
    with open(path) as f:
        data = json.load(f)

    # Domain
    domain_data = data.get("domain", {})
    domain = DomainConfig(
        width_um=domain_data.get("width_um", 1000),
        height_um=domain_data.get("height_um", 1000),
        env_grid_um=domain_data.get("env_grid_um", 10),
    )
    if domain.width_um % domain.env_grid_um != 0:
        raise ValueError(f"env_grid_um ({domain.env_grid_um}) must evenly divide width_um ({domain.width_um})")
    if domain.height_um % domain.env_grid_um != 0:
        raise ValueError(f"env_grid_um ({domain.env_grid_um}) must evenly divide height_um ({domain.height_um})")

    # Time
    time_data = data.get("time", {})
    time_cfg = TimeConfig(
        total_hours=time_data.get("total_hours", 168),
        dt_diffusion_hr=time_data.get("dt_diffusion_hr", 0.01),
        output_interval_hr=time_data.get("output_interval_hr", 1),
        burnin_hours=time_data.get("burnin_hours", 100),
    )

    # Environment
    env_data = data.get("environment", {})
    substrates = {}
    for name, sub_data in env_data.get("substrates", {}).items():
        substrates[name] = _parse_substrate(sub_data)
    environment = EnvironmentConfig(
        O2_blood_mmHg=env_data.get("O2_blood_mmHg", 60.0),
        q_O2_transfer_per_hr=env_data.get("q_O2_transfer_per_hr", 5.0),
        q_glucose_transfer_per_hr=env_data.get("q_glucose_transfer_per_hr", 3.0),
        hypoxia_threshold_mmHg=env_data.get("hypoxia_threshold_mmHg", 10.0),
        substrates=substrates,
        egfr_hypoxia_upregulation=env_data.get("egfr_hypoxia_upregulation", 1.5),
    )

    # Angiogenesis
    angio_data = data.get("angiogenesis", {})
    angiogenesis = AngiogenesisConfig(
        angiogenesis_threshold_nM=angio_data.get("angiogenesis_threshold_nM", 5.0),
        vessel_maturation_hr=angio_data.get("vessel_maturation_hr", 48.0),
        sprouting_rate_per_hr=angio_data.get("sprouting_rate_per_hr", 0.01),
        max_sprout_length=angio_data.get("max_sprout_length", 200),
    )

    # Immune recruitment
    immune_data = data.get("immune_recruitment", {})
    immune = ImmuneRecruitmentConfig(
        chemokine_threshold=immune_data.get("chemokine_threshold", 0.5),
        recruitment_rate_per_hr=immune_data.get("recruitment_rate_per_hr", 0.005),
        kill_radius_um=immune_data.get("kill_radius_um", 20.0),
        synapse_formation_time_hr=immune_data.get("synapse_formation_time_hr", 1.5),
        kill_probability_per_synapse=immune_data.get("kill_probability_per_synapse", 0.7),
        min_contact_for_kill_hr=immune_data.get("min_contact_for_kill_hr", 0.5),
        max_kills_before_exhaustion=immune_data.get("max_kills_before_exhaustion", 10),
        exhaustion_per_kill=immune_data.get("exhaustion_per_kill", 0.1),
        exhaustion_recovery_rate_per_hr=immune_data.get("exhaustion_recovery_rate_per_hr", 0.02),
        exhausted_kill_penalty=immune_data.get("exhausted_kill_penalty", 0.5),
        activation_radius_um=immune_data.get("activation_radius_um", 50.0),
        activation_rate_per_hr=immune_data.get("activation_rate_per_hr", 0.3),
        deactivation_rate_per_hr=immune_data.get("deactivation_rate_per_hr", 0.05),
        min_activation_for_kill=immune_data.get("min_activation_for_kill", 0.2),
        chemokine_diffusion_um2_hr=immune_data.get("chemokine_diffusion_um2_hr", 1800.0),
        chemokine_decay_per_hr=immune_data.get("chemokine_decay_per_hr", 0.5),
        tumor_chemokine_secretion=immune_data.get("tumor_chemokine_secretion", 1.0),
    )

    # Cell types
    cell_types: dict[int, CellTypeConfig] = {}
    for type_id_str, ct_data in data.get("cell_types", {}).items():
        type_id = int(type_id_str)
        cell_types[type_id] = _parse_cell_type(type_id, ct_data)

    # Tumor seeds
    tumor_seeds = _parse_tumor_seeds(data.get("tumor_seeds", []))

    # Initialization
    init_data = data.get("initialization", {})
    initialization = InitializationConfig(
        cell_density_per_um2=init_data.get("cell_density_per_um2", 0.005),
        fractions=init_data.get("fractions", {}),
        vascular_spacing_um=init_data.get("vascular_spacing_um", 150.0),
        astrocyte_domain_um=init_data.get("astrocyte_domain_um", 75.0),
    )

    return SimulationConfig(
        rng_seed=data.get("rng_seed", 42),
        output_dir=data.get("output_dir"),
        domain=domain,
        time=time_cfg,
        environment=environment,
        angiogenesis=angiogenesis,
        immune_recruitment=immune,
        cell_types=cell_types,
        tumor_seeds=tumor_seeds,
        initialization=initialization,
    )
