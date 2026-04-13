"""CAUSANTA Simulation Engine.

Generates synthetic tissue sections with known ground-truth causal structure
for benchmarking causal discovery algorithms.
"""

from .core import Simulation, main
from .config import SimulationConfig, load_config
from .cells import Cell, CellPopulation, create_cell
from .cells import NEURON, ASTROCYTE, OLIGODENDROCYTE, MICROGLIA
from .cells import ENDOTHELIAL, PERICYTE, TUMOR, RECRUITED_IMMUNE, NECROTIC

__all__ = [
    "Simulation",
    "main",
    "SimulationConfig",
    "load_config",
    "Cell",
    "CellPopulation",
    "create_cell",
    "NEURON",
    "ASTROCYTE",
    "OLIGODENDROCYTE",
    "MICROGLIA",
    "ENDOTHELIAL",
    "PERICYTE",
    "TUMOR",
    "RECRUITED_IMMUNE",
    "NECROTIC",
]
