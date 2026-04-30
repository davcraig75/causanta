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

# EGFR expression parameters (calibrated to match first-stage regression)
EGFR_BASE_EXPRESSION = 2.89   # Baseline EGFR expression (includes hypoxia effects in population)
EGFR_PER_ECDNA_COPY = 1.21    # Each ecDNA copy adds this much expression (gene dosage, κ)
EGFR_NOISE_CV = 0.1          # Coefficient of variation for transcriptional noise

# HIF-1alpha-driven EGFR upregulation under hypoxia. HIF binds the EGFR
# promoter and increases transcription. Reference: Franovic et al. 2007,
# Peng et al. 2006. This creates the confounder->exposure edge (M -> X)
# required for X to be endogenous w.r.t. hypoxia, which is what the SIV
# framework is designed to correct for. Without this term, corr(EGFR, M)=0
# by construction and IV has no OLS bias to correct.
EGFR_HYPOXIA_UPREGULATION = 1.5  # fractional boost under hypoxia (1+value)x mRNA


def compute_egfr_expression(
    ecDNA_count: int,
    is_hypoxic: bool = False,
    rng: np.random.Generator | None = None,
    include_noise: bool = True,
    hypoxia_upregulation: float | None = None,
) -> float:
    """Compute EGFR expression level from ecDNA copy number and hypoxia.

    The causal chain is:
        ecDNA_count (Z)  -> EGFR_mRNA -> EGFR_protein (X)
        hypoxia     (M)  -> HIF-2alpha -> EGFR_protein (additional endogenous driver)

    Gene dosage models additional EGFR transcripts from amplified copies on
    ecDNA. The hypoxia effect models translational upregulation of EGFR under
    hypoxia, as shown by Franovic et al. (2007) to be mediated by HIF-2alpha.

    Formula:
        EGFR = [base + kappa * ecDNA_count] * [1 + kappa_hyp * is_hypoxic] * noise

    Where noise ~ LogNormal(0, CV) captures transcriptional stochasticity.

    The hypoxia term is what makes EGFR endogenous w.r.t. the hypoxia
    confounder M. Without it, corr(EGFR, is_hypoxic) = 0 and IV has no
    omitted-variable bias to correct -- because the 'confounder' doesn't
    co-vary with the exposure. The strength of this effect can be varied
    for robustness analysis.

    Args:
        ecDNA_count: Number of ecDNA copies carrying EGFR
        is_hypoxic: Whether the cell is in a hypoxic microenvironment
        rng: Random number generator for noise (None = no noise)
        include_noise: Whether to add transcriptional noise
        hypoxia_upregulation: Fractional boost under hypoxia (default: 1.5 for 2.5x).
                              Set to 0 to remove hypoxia effect on EGFR.

    Returns:
        EGFR expression level (arbitrary units, ~1.0 for normal cells)

    Example (at ecDNA=20, hypoxia_upregulation=1.5):
        normoxic:  EGFR ~ 2.89 + 1.21*20 = 27.1
        hypoxic:   EGFR ~ (2.89 + 1.21*20) * 2.5 = 67.8

    References:
        Hung et al. (2021) Nature: ecDNA forms transcriptional hubs
        Franovic et al. (2007) PNAS: HIF-2alpha regulates EGFR translation
    """
    # Use default if not specified
    if hypoxia_upregulation is None:
        hypoxia_upregulation = EGFR_HYPOXIA_UPREGULATION

    # Cap ecDNA count
    capped_count = min(ecDNA_count, MAX_ECDNA_COPIES)

    # Gene dosage: more copies → more mRNA
    mean_expression = EGFR_BASE_EXPRESSION + EGFR_PER_ECDNA_COPY * capped_count

    # Hypoxia-driven upregulation (HIF-2alpha mediated translation)
    if is_hypoxic and hypoxia_upregulation > 0:
        mean_expression *= (1.0 + hypoxia_upregulation)

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


VEGF_NORMOXIC_BASAL_FRACTION = 0.2  # Normoxic cells secrete this fraction of their max rate


def modulate_vegf_secretion(
    base_rate: float,
    egfr_expression: float,
    beta: float,
    is_hypoxic: bool = False,
) -> float:
    """EGFR expression and hypoxia jointly drive VEGF secretion.

    Causal pathway:
        ecDNA (Z)  -> EGFR (X) -> VEGF transcription
        hypoxia(M) -> HIF-1alpha -> VEGF transcription (5-fold upregulation)

    EGFR effect is sqrt-saturating (receptor saturation). Hypoxia acts as
    a multiplicative switch: under normoxia cells secrete at a basal fraction
    (0.2x) of their EGFR-driven rate; HIF-1alpha fully activates secretion
    under hypoxia. Net: hypoxic cells secrete 5x more VEGF than normoxic
    at the same EGFR level.

    This matches the paper's structural equation:
        S_eff = S_base * (1 + beta * sqrt(X)) * (0.2 + 0.8 * M)

    Note: This function now returns the hypoxia-gated rate. Previously the
    gating was applied later in environment.py; it now lives here so that
    stored cell.VEGF_secretion reflects the actual secretion rate, which
    matters for IV/OLS analysis downstream.

    Args:
        base_rate: Base VEGF secretion rate (amol/hr)
        egfr_expression: EGFR expression level
        beta: EGFR effect size parameter
        is_hypoxic: Whether the cell is in a hypoxic microenvironment

    Returns:
        Effective VEGF secretion rate (amol/hr)
    """
    if egfr_expression <= 0:
        egfr_factor = 1.0
    else:
        egfr_factor = 1.0 + beta * np.sqrt(egfr_expression)
    hypoxia_factor = 1.0 if is_hypoxic else VEGF_NORMOXIC_BASAL_FRACTION
    return base_rate * egfr_factor * hypoxia_factor


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

    With gamma=0.5:
        - EGFR=1 (normal): a_eff = a_base / 1.5 (33% reduction)
        - EGFR=6 (10 ecDNA): a_eff = a_base / 2.4 (58% reduction)
        - EGFR=11 (20 ecDNA): a_eff = a_base / 2.8 (64% reduction)

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
