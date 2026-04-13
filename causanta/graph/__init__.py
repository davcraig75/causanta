"""Causal graph representation for CAUSANTA.

Provides explicit DAG data structures for representing ground-truth
causal models and comparing against discovered graphs.
"""

from .causal_dag import (
    CausalDAG,
    CausalEdge,
    CausalNode,
    NodeType,
    EffectType,
    build_causanta_dag,
    compare_dags,
)

__all__ = [
    "CausalDAG",
    "CausalEdge",
    "CausalNode",
    "NodeType",
    "EffectType",
    "build_causanta_dag",
    "compare_dags",
]
