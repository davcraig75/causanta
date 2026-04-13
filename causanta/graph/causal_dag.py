"""Explicit causal DAG representation for CAUSANTA.

Provides a formal graph structure for the ground-truth causal model,
supporting visualization, export, and comparison with discovered graphs.

This module is designed to be independent of the simulation engine,
accepting plain dictionaries for configuration to avoid circular imports.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any


class NodeType(Enum):
    """Type of node in the causal graph."""
    INSTRUMENT = "instrument"      # Randomized source (ecDNA segregation)
    EXPOSURE = "exposure"          # Instrumented variable (gene expression)
    MEDIATOR = "mediator"          # Intermediate variable
    OUTCOME = "outcome"            # Target phenotype
    CONFOUNDER = "confounder"      # Shared cause


class EffectType(Enum):
    """Functional form of causal effect."""
    LINEAR = "linear"              # Y = a + b*X
    LOG = "log"                    # Y = a / (1 + b*log2(1+X))
    MULTIPLICATIVE = "multiplicative"  # Y = a * (1 + b*X)
    INVERSE = "inverse"            # Y = a / (1 + b*X)


@dataclass
class CausalNode:
    """A node in the causal DAG."""

    name: str
    description: str
    node_type: NodeType
    units: str = ""
    observable: bool = True

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "description": self.description,
            "node_type": self.node_type.value,
            "units": self.units,
            "observable": self.observable,
        }


@dataclass
class CausalEdge:
    """A directed edge in the causal DAG representing a causal effect."""

    source: str
    target: str
    effect_type: EffectType
    parameter_symbol: str          # Greek letter (α, β, γ, δ)
    parameter_name: str            # Config key name
    parameter_value: float         # Actual value from config
    equation_template: str         # e.g., "T_eff = T_base / (1 + {param} * log2(1 + {source}))"
    description: str = ""

    @property
    def equation(self) -> str:
        """Format the equation with actual parameter value."""
        return self.equation_template.format(
            param=f"{self.parameter_symbol}={self.parameter_value:.3f}",
            source=self.source,
            target=self.target,
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "source": self.source,
            "target": self.target,
            "effect_type": self.effect_type.value,
            "parameter_symbol": self.parameter_symbol,
            "parameter_name": self.parameter_name,
            "parameter_value": self.parameter_value,
            "equation": self.equation,
            "description": self.description,
        }


@dataclass
class CausalDAG:
    """Complete causal directed acyclic graph for CAUSANTA simulation.

    The DAG encodes the ground-truth causal structure that the simulation
    implements. This allows:
    - Visualization of the causal model
    - Comparison with causal discovery algorithm outputs
    - Documentation of the true data-generating process
    """

    name: str
    description: str
    nodes: list[CausalNode] = field(default_factory=list)
    edges: list[CausalEdge] = field(default_factory=list)

    def add_node(self, node: CausalNode) -> None:
        """Add a node to the DAG."""
        self.nodes.append(node)

    def add_edge(self, edge: CausalEdge) -> None:
        """Add an edge to the DAG."""
        self.edges.append(edge)

    def get_node(self, name: str) -> CausalNode | None:
        """Get a node by name."""
        for node in self.nodes:
            if node.name == name:
                return node
        return None

    def get_parents(self, node_name: str) -> list[str]:
        """Get the names of all parent nodes."""
        return [e.source for e in self.edges if e.target == node_name]

    def get_children(self, node_name: str) -> list[str]:
        """Get the names of all child nodes."""
        return [e.target for e in self.edges if e.source == node_name]

    def get_edge(self, source: str, target: str) -> CausalEdge | None:
        """Get an edge by source and target."""
        for edge in self.edges:
            if edge.source == source and edge.target == target:
                return edge
        return None

    def to_dict(self) -> dict[str, Any]:
        """Convert DAG to dictionary for JSON serialization."""
        return {
            "name": self.name,
            "description": self.description,
            "nodes": [n.to_dict() for n in self.nodes],
            "edges": [e.to_dict() for e in self.edges],
        }

    def to_mermaid(self, direction: str = "TB") -> str:
        """Generate Mermaid diagram syntax for visualization.

        Args:
            direction: Graph direction (TB=top-bottom, LR=left-right)

        Returns:
            Mermaid diagram as string.
        """
        lines = [f"graph {direction}"]

        # Define node styles by type
        style_map = {
            NodeType.INSTRUMENT: "fill:#e8f5e9,stroke:#2e7d32",
            NodeType.EXPOSURE: "fill:#e3f2fd,stroke:#1565c0",
            NodeType.MEDIATOR: "fill:#fff3e0,stroke:#e65100",
            NodeType.OUTCOME: "fill:#fce4ec,stroke:#c2185b",
            NodeType.CONFOUNDER: "fill:#f3e5f5,stroke:#7b1fa2",
        }

        # Add nodes with styling
        for node in self.nodes:
            node_id = node.name.replace(" ", "_")
            lines.append(f"    {node_id}[{node.name}]")

        # Add edges with labels
        for edge in self.edges:
            src = edge.source.replace(" ", "_")
            tgt = edge.target.replace(" ", "_")
            label = f"{edge.parameter_symbol}={edge.parameter_value:.2f}"
            lines.append(f"    {src} -->|{label}| {tgt}")

        # Add style classes
        lines.append("")
        for node in self.nodes:
            node_id = node.name.replace(" ", "_")
            style = style_map.get(node.node_type, "fill:#fff")
            lines.append(f"    style {node_id} {style}")

        return "\n".join(lines)

    def to_graphviz(self) -> str:
        """Generate Graphviz DOT syntax for visualization."""
        lines = ["digraph CausalDAG {"]
        lines.append("    rankdir=TB;")
        lines.append("    node [shape=box, style=filled];")

        # Color map by node type
        color_map = {
            NodeType.INSTRUMENT: "lightgreen",
            NodeType.EXPOSURE: "lightblue",
            NodeType.MEDIATOR: "lightyellow",
            NodeType.OUTCOME: "lightpink",
            NodeType.CONFOUNDER: "lavender",
        }

        # Add nodes
        for node in self.nodes:
            color = color_map.get(node.node_type, "white")
            label = f"{node.name}\\n({node.units})" if node.units else node.name
            lines.append(f'    "{node.name}" [label="{label}", fillcolor={color}];')

        # Add edges
        for edge in self.edges:
            label = f"{edge.parameter_symbol}={edge.parameter_value:.2f}"
            lines.append(f'    "{edge.source}" -> "{edge.target}" [label="{label}"];')

        lines.append("}")
        return "\n".join(lines)

    def to_adjacency_matrix(self) -> tuple[list[str], list[list[float]]]:
        """Convert to adjacency matrix with edge weights.

        Returns:
            (node_names, adjacency_matrix) where matrix[i][j] is the
            parameter value for edge i->j (0 if no edge).
        """
        node_names = [n.name for n in self.nodes]
        n = len(node_names)
        name_to_idx = {name: i for i, name in enumerate(node_names)}

        matrix = [[0.0] * n for _ in range(n)]

        for edge in self.edges:
            if edge.source in name_to_idx and edge.target in name_to_idx:
                i = name_to_idx[edge.source]
                j = name_to_idx[edge.target]
                matrix[i][j] = edge.parameter_value

        return node_names, matrix


def build_causanta_dag(params: dict[str, Any] | None = None) -> CausalDAG:
    """Build the CAUSANTA causal DAG from parameter dictionary.

    This constructs the ground-truth causal graph that the simulation
    implements, extracting effect sizes from the provided parameters.

    Args:
        params: Dictionary with effect parameters. If None, uses defaults.
            Expected keys:
            - alpha (ecDNA_effect_on_division): float
            - beta (ecDNA_effect_on_VEGF): float
            - delta (ecDNA_effect_on_migration): float
            - gamma (ecDNA_effect_on_survival): float
            - angiogenesis_threshold_nM: float
            - q_O2_transfer_per_hr: float
            - O2_prolif_threshold_mmHg: float
            - glucose_prolif_threshold_mM: float

    Returns:
        CausalDAG with complete CAUSANTA model structure.
    """
    if params is None:
        params = {}

    # Default parameter values (match simulation defaults)
    alpha = params.get("alpha", params.get("ecDNA_effect_on_division", 0.5))
    beta = params.get("beta", params.get("ecDNA_effect_on_VEGF", 0.02))
    delta = params.get("delta", params.get("ecDNA_effect_on_migration", 0.01))
    gamma = params.get("gamma", params.get("ecDNA_effect_on_survival", 0.1))
    angio_threshold = params.get("angiogenesis_threshold_nM", 1.0)
    q_o2 = params.get("q_O2_transfer_per_hr", 0.0005)
    o2_threshold = params.get("O2_prolif_threshold_mmHg", 1.0)
    glucose_threshold = params.get("glucose_prolif_threshold_mM", 0.5)

    dag = CausalDAG(
        name="CAUSANTA Glioma Model",
        description=(
            "Ground-truth causal structure for glioma tumor simulation. "
            "ecDNA copy number acts as a somatic instrumental variable (SIV), "
            "with binomial segregation during mitosis providing natural randomization."
        ),
    )

    # === NODES ===

    # Instrument (randomized by mitosis)
    dag.add_node(CausalNode(
        name="ecDNA_count",
        description="Extrachromosomal DNA copy number per cell",
        node_type=NodeType.INSTRUMENT,
        units="copies",
        observable=True,
    ))

    # Exposures (directly affected by instrument)
    dag.add_node(CausalNode(
        name="EGFR_expression",
        description="EGFR gene expression level (carried on ecDNA)",
        node_type=NodeType.EXPOSURE,
        units="relative",
        observable=True,
    ))

    dag.add_node(CausalNode(
        name="MYC_expression",
        description="MYC gene expression level (carried on ecDNA)",
        node_type=NodeType.EXPOSURE,
        units="relative",
        observable=True,
    ))

    # Mediators
    dag.add_node(CausalNode(
        name="VEGF_secretion",
        description="Vascular endothelial growth factor secretion rate",
        node_type=NodeType.MEDIATOR,
        units="amol/hr",
        observable=True,
    ))

    dag.add_node(CausalNode(
        name="angiogenesis",
        description="Local vascular sprouting and vessel density",
        node_type=NodeType.MEDIATOR,
        units="density",
        observable=True,
    ))

    # Outcomes
    dag.add_node(CausalNode(
        name="proliferation_rate",
        description="Cell division rate (inverse of cycle time)",
        node_type=NodeType.OUTCOME,
        units="1/hr",
        observable=True,
    ))

    dag.add_node(CausalNode(
        name="migration_speed",
        description="Cell migration velocity",
        node_type=NodeType.OUTCOME,
        units="um/hr",
        observable=True,
    ))

    dag.add_node(CausalNode(
        name="apoptosis_resistance",
        description="Resistance to programmed cell death",
        node_type=NodeType.OUTCOME,
        units="relative",
        observable=False,  # Inferred from survival
    ))

    # Confounders
    dag.add_node(CausalNode(
        name="O2_local",
        description="Local oxygen partial pressure",
        node_type=NodeType.CONFOUNDER,
        units="mmHg",
        observable=True,
    ))

    dag.add_node(CausalNode(
        name="glucose_local",
        description="Local glucose concentration",
        node_type=NodeType.CONFOUNDER,
        units="mM",
        observable=True,
    ))

    # === EDGES ===

    # ecDNA -> gene expression (implicit, always active)
    dag.add_edge(CausalEdge(
        source="ecDNA_count",
        target="EGFR_expression",
        effect_type=EffectType.LINEAR,
        parameter_symbol="κ",
        parameter_name="ecDNA_cargo",
        parameter_value=1.0,  # Direct relationship
        equation_template="EGFR = {param} * ecDNA_count (cargo gene)",
        description="ecDNA carries EGFR amplicon; copy number directly determines expression",
    ))

    dag.add_edge(CausalEdge(
        source="ecDNA_count",
        target="MYC_expression",
        effect_type=EffectType.LINEAR,
        parameter_symbol="κ",
        parameter_name="ecDNA_cargo",
        parameter_value=1.0,
        equation_template="MYC = {param} * ecDNA_count (cargo gene)",
        description="ecDNA carries MYC amplicon; copy number directly determines expression",
    ))

    # ecDNA/expression -> proliferation (via α)
    dag.add_edge(CausalEdge(
        source="ecDNA_count",
        target="proliferation_rate",
        effect_type=EffectType.LOG,
        parameter_symbol="α",
        parameter_name="ecDNA_effect_on_division",
        parameter_value=alpha,
        equation_template="T_eff = T_base / (1 + {param} * log2(1 + ecDNA))",
        description="Higher ecDNA accelerates cell cycle via oncogene overexpression",
    ))

    # ecDNA -> VEGF secretion (via β)
    dag.add_edge(CausalEdge(
        source="ecDNA_count",
        target="VEGF_secretion",
        effect_type=EffectType.MULTIPLICATIVE,
        parameter_symbol="β",
        parameter_name="ecDNA_effect_on_VEGF",
        parameter_value=beta,
        equation_template="S_eff = S_base * (1 + {param} * ecDNA)",
        description="Higher ecDNA increases VEGF production under hypoxia",
    ))

    # ecDNA -> migration (via δ)
    dag.add_edge(CausalEdge(
        source="ecDNA_count",
        target="migration_speed",
        effect_type=EffectType.MULTIPLICATIVE,
        parameter_symbol="δ",
        parameter_name="ecDNA_effect_on_migration",
        parameter_value=delta,
        equation_template="v_eff = v_base * (1 + {param} * ecDNA)",
        description="Higher ecDNA increases invasive migration capacity",
    ))

    # ecDNA -> apoptosis resistance (via γ)
    dag.add_edge(CausalEdge(
        source="ecDNA_count",
        target="apoptosis_resistance",
        effect_type=EffectType.INVERSE,
        parameter_symbol="γ",
        parameter_name="ecDNA_effect_on_survival",
        parameter_value=gamma,
        equation_template="a_eff = a_base / (1 + {param} * ecDNA)",
        description="Higher ecDNA confers resistance to apoptosis",
    ))

    # VEGF -> angiogenesis
    dag.add_edge(CausalEdge(
        source="VEGF_secretion",
        target="angiogenesis",
        effect_type=EffectType.LINEAR,
        parameter_symbol="θ",
        parameter_name="angiogenesis_threshold_nM",
        parameter_value=angio_threshold,
        equation_template="sprouting when VEGF > {param}",
        description="VEGF above threshold triggers endothelial sprouting",
    ))

    # angiogenesis -> O2 (new vessels increase oxygen)
    dag.add_edge(CausalEdge(
        source="angiogenesis",
        target="O2_local",
        effect_type=EffectType.LINEAR,
        parameter_symbol="q",
        parameter_name="q_O2_transfer_per_hr",
        parameter_value=q_o2,
        equation_template="O2_supply = {param} * vascular_density * (O2_blood - O2_local)",
        description="New vessels increase local oxygen delivery",
    ))

    # O2 -> proliferation (confounder path)
    dag.add_edge(CausalEdge(
        source="O2_local",
        target="proliferation_rate",
        effect_type=EffectType.LINEAR,
        parameter_symbol="τ",
        parameter_name="O2_prolif_threshold_mmHg",
        parameter_value=o2_threshold,
        equation_template="proliferation requires O2 > {param}",
        description="Oxygen is required for cell cycle progression",
    ))

    # glucose -> proliferation (confounder path)
    dag.add_edge(CausalEdge(
        source="glucose_local",
        target="proliferation_rate",
        effect_type=EffectType.LINEAR,
        parameter_symbol="ρ",
        parameter_name="glucose_prolif_threshold_mM",
        parameter_value=glucose_threshold,
        equation_template="proliferation requires glucose > {param}",
        description="Glucose is required for cell cycle progression",
    ))

    return dag


def compare_dags(ground_truth: CausalDAG, discovered: dict[str, list[str]]) -> dict[str, Any]:
    """Compare a discovered causal graph against ground truth.

    Args:
        ground_truth: The true causal DAG from the simulation
        discovered: Dictionary mapping each node to its discovered parents

    Returns:
        Metrics including true positives, false positives, false negatives,
        precision, recall, F1 score, and structural Hamming distance.
    """
    # Build ground truth edge set
    gt_edges = {(e.source, e.target) for e in ground_truth.edges}

    # Build discovered edge set
    disc_edges = set()
    for target, parents in discovered.items():
        for parent in parents:
            disc_edges.add((parent, target))

    # Compute metrics
    true_positives = gt_edges & disc_edges
    false_positives = disc_edges - gt_edges
    false_negatives = gt_edges - disc_edges

    tp = len(true_positives)
    fp = len(false_positives)
    fn = len(false_negatives)

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0

    # Structural Hamming Distance (SHD)
    shd = fp + fn

    return {
        "true_positives": list(true_positives),
        "false_positives": list(false_positives),
        "false_negatives": list(false_negatives),
        "n_true_positives": tp,
        "n_false_positives": fp,
        "n_false_negatives": fn,
        "precision": precision,
        "recall": recall,
        "f1_score": f1,
        "structural_hamming_distance": shd,
        "n_gt_edges": len(gt_edges),
        "n_discovered_edges": len(disc_edges),
    }
