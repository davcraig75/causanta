"""ecDNA segregation and phenotype modulation for CAUSANTA.

Implements the Somatic Instrumental Variable (SIV) mechanism:
binomial segregation during mitosis and ground-truth causal effects
of ecDNA copy number on cell phenotypes.
"""

from __future__ import annotations

import numpy as np


def segregate_ecdna(
    parent_count: int,
    p: float,
    rng: np.random.Generator,
) -> tuple[int, int]:
    """Binomial segregation of ecDNA during mitosis.

    Each ecDNA copy independently goes to the daughter with probability p.

    Returns:
        (parent_new_count, daughter_count)
    """
    if parent_count <= 0:
        return 0, 0
    daughter = int(rng.binomial(parent_count, p))
    return parent_count - daughter, daughter


def modulate_division_time(
    base_time_hr: float,
    ecDNA_count: int,
    alpha: float,
) -> float:
    """ecDNA accelerates proliferation.

    T_eff = T_base / (1 + alpha * log2(1 + ecDNA_count))
    """
    if alpha == 0.0 or ecDNA_count == 0:
        return base_time_hr
    return base_time_hr / (1.0 + alpha * np.log2(1.0 + ecDNA_count))


def modulate_vegf_secretion(
    base_rate: float,
    ecDNA_count: int,
    beta: float,
) -> float:
    """ecDNA amplifies VEGF secretion.

    S_eff = S_base * (1 + beta * ecDNA_count)
    """
    return base_rate * (1.0 + beta * ecDNA_count)


def modulate_migration_speed(
    base_speed: float,
    ecDNA_count: int,
    delta: float,
) -> float:
    """ecDNA increases migration speed.

    v_eff = v_base * (1 + delta * ecDNA_count)
    """
    return base_speed * (1.0 + delta * ecDNA_count)


def modulate_apoptosis_rate(
    base_rate: float,
    ecDNA_count: int,
    gamma: float,
) -> float:
    """ecDNA confers apoptosis resistance.

    a_eff = a_base / (1 + gamma * ecDNA_count)
    """
    if gamma == 0.0 or ecDNA_count == 0:
        return base_rate
    return base_rate / (1.0 + gamma * ecDNA_count)
