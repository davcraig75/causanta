"""ecDNA segregation and phenotype modulation for CAUSANTA.

Implements the Somatic Instrumental Variable (SIV) mechanism:
ecDNA replication during S-phase followed by binomial segregation
during mitosis, with ground-truth causal effects of ecDNA copy
number on cell phenotypes.

References:
- Nature 2025: "Genetic elements promote retention of extrachromosomal DNA"
- Nature Genetics 2022: "Evolutionary dynamics of extrachromosomal DNA"
- Cell 2020: "ecDNA promotes tumor heterogeneity through random inheritance"

Key biology:
- ecDNA lacks centromeres → random segregation during mitosis
- ecDNA replicates during S-phase → approximately doubles before division
- EGFR amplification on ecDNA → increased proliferation and invasion
"""

from __future__ import annotations

import numpy as np

# Maximum ecDNA copy number (biological constraint - very high copy numbers
# are associated with genomic instability and reduced fitness)
MAX_ECDNA_COPIES = 100

# EGFR expression parameters
EGFR_BASE_EXPRESSION = 1.0  # Normal cells have baseline EGFR expression
EGFR_PER_ECDNA_COPY = 0.5   # Each ecDNA copy adds this much expression
EGFR_NOISE_CV = 0.1         # Coefficient of variation for transcriptional noise


def compute_egfr_expression(
    ecDNA_count: int,
    rng: np.random.Generator | None = None,
    include_noise: bool = True,
) -> float:
    """Compute EGFR expression level from ecDNA copy number.

    The causal chain is:
        ecDNA_count → EGFR_mRNA → EGFR_protein (expression)

    This function models gene dosage effects: more ecDNA copies carrying
    EGFR → more EGFR transcription → higher EGFR protein levels.

    Formula:
        EGFR = base + (copies_per_ecDNA * ecDNA_count) * noise

    Where noise ~ LogNormal(0, CV) to model transcriptional stochasticity.

    Args:
        ecDNA_count: Number of ecDNA copies carrying EGFR
        rng: Random number generator for noise (None = no noise)
        include_noise: Whether to add transcriptional noise

    Returns:
        EGFR expression level (arbitrary units, ~1.0 for normal cells)

    Example:
        ecDNA=0:  EGFR ≈ 1.0 (baseline)
        ecDNA=10: EGFR ≈ 1.0 + 0.5*10 = 6.0 (6x overexpression)
        ecDNA=20: EGFR ≈ 1.0 + 0.5*20 = 11.0 (11x overexpression)
        ecDNA=50: EGFR ≈ 1.0 + 0.5*50 = 26.0 (26x overexpression)

    References:
        Hung et al. (2021) Nature: ecDNA forms transcriptional hubs
    """
    # Cap ecDNA count
    capped_count = min(ecDNA_count, MAX_ECDNA_COPIES)

    # Gene dosage: more copies → more expression
    mean_expression = EGFR_BASE_EXPRESSION + EGFR_PER_ECDNA_COPY * capped_count

    # Add transcriptional noise if requested
    if include_noise and rng is not None and EGFR_NOISE_CV > 0:
        # Log-normal noise preserves positivity and gives realistic variance
        noise = rng.lognormal(0, EGFR_NOISE_CV)
        expression = mean_expression * noise
    else:
        expression = mean_expression

    return max(0.0, expression)  # Expression can't be negative


def replicate_and_segregate_ecdna(
    parent_count: int,
    p: float,
    rng: np.random.Generator,
    replication_fidelity: float = 0.95,
    max_copies: int = MAX_ECDNA_COPIES,
) -> tuple[int, int]:
    """ecDNA replication during S-phase followed by binomial segregation.

    Biological model:
    1. S-phase: ecDNA replicates (approximately doubles, with some noise)
    2. M-phase: Replicated ecDNAs segregate randomly to daughters

    The replication step ensures ecDNA copy number is maintained across
    generations (in expectation), while segregation creates variance that
    enables natural selection to act on ecDNA-driven phenotypes.

    Args:
        parent_count: Number of ecDNA copies before replication
        p: Segregation probability (probability each copy goes to daughter)
        rng: Random number generator
        replication_fidelity: Probability each ecDNA successfully replicates
                             (default 0.95, accounts for replication errors)
        max_copies: Maximum allowed ecDNA copies (caps runaway amplification)

    Returns:
        (parent_new_count, daughter_count) after segregation
    """
    if parent_count <= 0:
        return 0, 0

    # Step 1: Replication during S-phase
    # Each ecDNA copy replicates with probability = replication_fidelity
    # This gives ~2x copies on average (with some variance)
    replicated_count = parent_count  # Original copies
    for _ in range(parent_count):
        if rng.random() < replication_fidelity:
            replicated_count += 1

    # Cap at maximum (genomic instability limits)
    replicated_count = min(replicated_count, max_copies)

    # Step 2: Random segregation during M-phase
    # Each of the replicated copies goes to daughter with probability p
    daughter = int(rng.binomial(replicated_count, p))
    parent_new = replicated_count - daughter

    # Apply maximum cap to both
    parent_new = min(parent_new, max_copies)
    daughter = min(daughter, max_copies)

    return parent_new, daughter


def segregate_ecdna(
    parent_count: int,
    p: float,
    rng: np.random.Generator,
) -> tuple[int, int]:
    """ecDNA replication and segregation during cell division.

    This is the main entry point that models the complete ecDNA
    inheritance process: replication followed by random segregation.

    Returns:
        (parent_new_count, daughter_count)
    """
    return replicate_and_segregate_ecdna(parent_count, p, rng)


def modulate_division_time(
    base_time_hr: float,
    egfr_expression: float,
    alpha: float,
) -> float:
    """EGFR expression accelerates proliferation.

    Causal pathway:
        ecDNA_count → EGFR_expression → RAS/MAPK, PI3K/AKT → cell cycle

    EGFR signaling promotes cell cycle progression. Higher EGFR expression
    → stronger mitogenic signaling → faster division.

    Formula (saturating effect):
        T_eff = T_base / (1 + alpha * log2(1 + EGFR_expression))

    With alpha=0.3:
        - EGFR=1 (normal): T_eff = T_base / 1.3 = 27.7h
        - EGFR=6 (10 ecDNA): T_eff = 36 / 2.15 = 16.7h
        - EGFR=11 (20 ecDNA): T_eff = 36 / 2.47 = 14.6h
        - EGFR=26 (50 ecDNA): T_eff = 36 / 2.87 = 12.5h

    Args:
        base_time_hr: Base division time for this cell type
        egfr_expression: EGFR expression level (from compute_egfr_expression)
        alpha: Effect size parameter

    Returns:
        Effective division time in hours
    """
    if alpha == 0.0 or egfr_expression <= 0:
        return base_time_hr
    return base_time_hr / (1.0 + alpha * np.log2(1.0 + egfr_expression))


def modulate_vegf_secretion(
    base_rate: float,
    egfr_expression: float,
    beta: float,
) -> float:
    """EGFR expression increases VEGF secretion capacity.

    Causal pathway:
        ecDNA → EGFR_expression → HIF-1α stabilization → VEGF transcription

    EGFR signaling can enhance HIF-1α stability and VEGF transcription,
    promoting angiogenesis. Effect has diminishing returns (sqrt).

    Formula:
        S_eff = S_base * (1 + beta * sqrt(EGFR_expression))

    With beta=0.1 and base=600:
        - EGFR=1 (normal): S_eff = 600 * 1.1 = 660
        - EGFR=6 (10 ecDNA): S_eff = 600 * 1.24 = 747
        - EGFR=11 (20 ecDNA): S_eff = 600 * 1.33 = 798
        - EGFR=26 (50 ecDNA): S_eff = 600 * 1.51 = 906

    Args:
        base_rate: Base VEGF secretion rate
        egfr_expression: EGFR expression level
        beta: Effect size parameter

    Returns:
        Effective VEGF secretion rate (amol/hr)
    """
    if egfr_expression <= 0:
        return base_rate
    return base_rate * (1.0 + beta * np.sqrt(egfr_expression))


def modulate_migration_speed(
    base_speed: float,
    egfr_expression: float,
    delta: float,
    is_hypoxic: bool = False,
    hypoxia_invasion_boost: float = 2.0,
) -> float:
    """Migration speed modulated by EGFR expression AND hypoxia ("Go or Grow").

    Causal pathway:
        ecDNA → EGFR_expression → EMT, MMPs → invasion capacity
        Hypoxia → HIF-1α → additional invasion boost

    Implements the "Go or Grow" hypothesis from glioblastoma biology:
    - Under normoxia: EGFR drives proliferation, modest migration
    - Under hypoxia: Cells switch to invasive phenotype ("escape")

    Two effects combine multiplicatively:
    1. EGFR expression increases baseline invasion capacity (EMT, MMPs)
    2. Hypoxia triggers "Go" response via HIF-1α pathway

    Formula:
        v_eff = v_base * (1 + delta * EGFR_expression) * hypoxia_multiplier

    With delta=0.05, base=10 um/hr, hypoxia_boost=2.0:
        Normoxic, EGFR=1:  10 * 1.05 * 1.0 = 10.5 um/hr
        Normoxic, EGFR=11: 10 * 1.55 * 1.0 = 15.5 um/hr
        Hypoxic,  EGFR=11: 10 * 1.55 * 2.0 = 31.0 um/hr

    Args:
        base_speed: Base migration speed (um/hr)
        egfr_expression: EGFR expression level
        delta: Effect size (invasion coefficient)
        is_hypoxic: Whether cell is in hypoxic environment
        hypoxia_invasion_boost: Fold-increase in migration under hypoxia

    Returns:
        Effective migration speed (um/hr)

    References:
        Kathagen-Buhmann et al. (2016) Neuro-Oncology: Go-or-Grow hypothesis
        Li et al. (2022) Nature Cell Biology: EGFR context-dependent effects
    """
    egfr_effect = 1.0 + delta * egfr_expression
    hypoxia_multiplier = hypoxia_invasion_boost if is_hypoxic else 1.0
    return base_speed * egfr_effect * hypoxia_multiplier


def modulate_apoptosis_rate(
    base_rate: float,
    egfr_expression: float,
    gamma: float,
) -> float:
    """EGFR expression confers apoptosis resistance.

    Causal pathway:
        ecDNA → EGFR_expression → PI3K/AKT → BAD/FOXO inhibition → survival

    EGFR signaling activates PI3K/AKT pathway, which phosphorylates and
    inhibits pro-apoptotic proteins, promoting cell survival.

    Formula (saturating protection):
        a_eff = a_base / (1 + gamma * log2(1 + EGFR_expression))

    With gamma=0.2:
        - EGFR=1 (normal): a_eff = a_base / 1.2 (17% reduction)
        - EGFR=6 (10 ecDNA): a_eff = a_base / 1.56 (36% reduction)
        - EGFR=11 (20 ecDNA): a_eff = a_base / 1.72 (42% reduction)

    Args:
        base_rate: Base apoptosis rate
        egfr_expression: EGFR expression level
        gamma: Effect size parameter

    Returns:
        Effective apoptosis rate (per hour)
    """
    if gamma == 0.0 or egfr_expression <= 0:
        return base_rate
    return base_rate / (1.0 + gamma * np.log2(1.0 + egfr_expression))
