"""CAUSANTA: Causal Analysis Using Somatic And Neighborhood Tissue Architecture.

CAUSANTA is a framework for:
1. Simulating synthetic tissue with known causal structure (causanta.simulate)
2. Analyzing causal effects from spatial data (causanta.analyze)

The simulation generates ground-truth data where ecDNA acts as a somatic
instrumental variable, enabling rigorous benchmarking of causal discovery algorithms.
"""

__version__ = "0.2.0"

# Import subpackages for convenient access
from . import simulate
from . import analyze
from . import graph

# Re-export key classes at top level for convenience
from .graph import CausalDAG, CausalEdge, CausalNode, build_causanta_dag
from .simulate import Simulation, SimulationConfig

__all__ = [
    # Subpackages
    "simulate",
    "analyze",
    "graph",
    # Key classes
    "CausalDAG",
    "CausalEdge",
    "CausalNode",
    "build_causanta_dag",
    "Simulation",
    "SimulationConfig",
]
