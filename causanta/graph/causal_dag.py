"""Explicit causal DAG representation for CAUSANTA.

Provides a formal graph structure for the ground-truth causal model,
supporting visualization, export, and comparison with discovered graphs.

This module is designed to be independent of the simulation engine,
accepting plain dictionaries for configuration to avoid circular imports.

Literature-supported causal structure:
- Lange et al. (2022) Nature Genetics: ecDNA segregates binomially
- Li et al. (2022) Nature Cell Biology: EGFR ligand context affects phenotype
- Kathagen-Buhmann et al. (2016) Neuro-Oncology: Hypoxia drives go-or-grow
- Hung et al. (2021) Nature: ecDNA forms transcriptional hubs
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
    SQRT = "sqrt"                  # Y = a * (1 + b*sqrt(X))
    MULTIPLICATIVE = "multiplicative"  # Y = a * (1 + b*X)
    INVERSE = "inverse"            # Y = a / (1 + b*X)
    THRESHOLD = "threshold"        # Y = f(X) if X > threshold


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
    sign: str = "+"                # Direction of effect: "+" or "-"

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
            "sign": self.sign,
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
            sign = edge.sign
            label = f"{sign}{edge.parameter_symbol}"
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
            label = f"{edge.sign}{edge.parameter_symbol}"
            color = "red" if edge.sign == "-" else "black"
            lines.append(f'    "{edge.source}" -> "{edge.target}" [label="{label}", color={color}];')

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
                # Use signed value
                sign = -1.0 if edge.sign == "-" else 1.0
                matrix[i][j] = sign * edge.parameter_value

        return node_names, matrix


def build_causanta_dag(params: dict[str, Any] | None = None) -> CausalDAG:
    """Build the CAUSANTA causal DAG from parameter dictionary.

    This constructs the ground-truth causal graph based on literature-supported
    relationships in glioblastoma:

    Literature basis:
    - Lange et al. (2022) Nature Genetics: ecDNA binomial segregation
    - Li et al. (2022) Nature Cell Biology: EGFR context-dependent effects
    - Kathagen-Buhmann et al. (2016) Neuro-Oncology: Hypoxia go-or-grow
    - Hung et al. (2021) Nature: ecDNA transcriptional hubs

    Core causal structure ("Go or Grow" model):
    1. ecDNA (EGFR) → Proliferation (+): EGFR signaling drives cell cycle
    2. Proliferation → O2 consumption: Growing cells deplete oxygen
    3. Low O2 = Hypoxia: Creates metabolic stress
    4. Hypoxia → Invasion (+): "Go" response - cells escape hypoxic core
    5. Hypoxia → Proliferation (-): "Grow" suppressed under hypoxia
    6. Hypoxia → VEGF (+): HIF-1α induces angiogenic signaling
    7. VEGF → Vessels (+): Angiogenesis
    8. Vessels → O2 (+): Blood delivers oxygen

    Args:
        params: Dictionary with effect parameters. If None, uses defaults.

    Returns:
        CausalDAG with complete CAUSANTA model structure.
    """
    if params is None:
        params = {}

    # Extract parameters (match config keys)
    alpha = params.get("alpha", params.get("ecDNA_effect_on_division", 0.3))
    beta = params.get("beta", params.get("ecDNA_effect_on_VEGF", 0.1))
    delta = params.get("delta", params.get("ecDNA_effect_on_migration", 0.05))
    gamma = params.get("gamma", params.get("ecDNA_effect_on_survival", 0.2))

    # EGFR expression parameters (from ecdna.py)
    egfr_per_copy = params.get("EGFR_per_ecDNA_copy", 0.5)
    egfr_base = params.get("EGFR_base_expression", 1.0)

    # Microenvironment parameters
    o2_consumption = params.get("O2_consumption_amol_hr", 72000.0)
    hypoxia_threshold = params.get("hypoxia_threshold_mmHg", 10.0)
    angio_threshold = params.get("angiogenesis_threshold_nM", 5.0)
    q_o2 = params.get("q_O2_transfer_per_hr", 5.0)

    dag = CausalDAG(
        name="CAUSANTA Glioblastoma Model",
        description=(
            "Ground-truth causal structure for glioblastoma simulation implementing "
            "the 'Go or Grow' hypothesis. ecDNA copy number (carrying EGFR) acts as "
            "a somatic instrumental variable (SIV) with binomial segregation providing "
            "natural randomization. Causal chain: ecDNA_count → EGFR_expression → Phenotypes. "
            "Key pathways: (1) ecDNA → EGFR → Proliferation, "
            "(2) Hypoxia → Invasion, (3) Proliferation → O2 depletion → Hypoxia feedback. "
            "Ref: Hung et al. 2021 (ecDNA transcriptional hubs)"
        ),
    )

    # =========================================================================
    # NODES
    # =========================================================================

    # --- Instrument (randomized by mitosis) ---
    dag.add_node(CausalNode(
        name="ecDNA_EGFR",
        description=(
            "Extrachromosomal DNA copy number carrying EGFR amplicon. "
            "Segregates binomially during mitosis (lacks centromere). "
            "Ref: Lange et al. 2022 Nature Genetics"
        ),
        node_type=NodeType.INSTRUMENT,
        units="copies (0-100)",
        observable=True,
    ))

    # --- Exposure (instrumented by ecDNA) ---
    dag.add_node(CausalNode(
        name="EGFR_expression",
        description=(
            "EGFR protein expression level, determined by ecDNA copy number "
            "(gene dosage effect) plus transcriptional noise. This is the "
            "exposure in the SIV framework - what we actually measure. "
            "Ref: Hung et al. 2021 Nature (ecDNA forms transcriptional hubs)"
        ),
        node_type=NodeType.EXPOSURE,
        units="arbitrary units",
        observable=True,
    ))

    # --- Primary Outcomes ---
    dag.add_node(CausalNode(
        name="Proliferation",
        description=(
            "Cell division rate (inverse of cycle time). Driven by EGFR "
            "signaling through RAS/MAPK and PI3K/AKT pathways. Suppressed "
            "under hypoxia due to metabolic constraints."
        ),
        node_type=NodeType.OUTCOME,
        units="1/hr",
        observable=True,
    ))

    dag.add_node(CausalNode(
        name="Invasion",
        description=(
            "Cell migration/invasion capacity. Activated under hypoxia as "
            "'escape' response (Go or Grow hypothesis). Involves EMT, MMP "
            "secretion, cytoskeletal reorganization. "
            "Ref: Kathagen-Buhmann et al. 2016"
        ),
        node_type=NodeType.OUTCOME,
        units="um/hr",
        observable=True,
    ))

    # --- Microenvironment (Mediators/Confounders) ---
    dag.add_node(CausalNode(
        name="O2_local",
        description=(
            "Local oxygen partial pressure. Depleted by proliferating cells, "
            "replenished by vasculature. Determines hypoxic state."
        ),
        node_type=NodeType.CONFOUNDER,
        units="mmHg",
        observable=True,
    ))

    dag.add_node(CausalNode(
        name="Hypoxia",
        description=(
            "Hypoxic state (O2 < threshold). Activates HIF-1α leading to "
            "metabolic reprogramming, VEGF secretion, and phenotype switch "
            "from proliferation to invasion."
        ),
        node_type=NodeType.MEDIATOR,
        units="binary/graded",
        observable=True,
    ))

    dag.add_node(CausalNode(
        name="VEGF",
        description=(
            "Vascular endothelial growth factor. Secreted under hypoxia "
            "(HIF-1α target gene). Promotes angiogenesis."
        ),
        node_type=NodeType.MEDIATOR,
        units="nM",
        observable=True,
    ))

    dag.add_node(CausalNode(
        name="Vasculature",
        description=(
            "Local blood vessel density. Source of oxygen and nutrients. "
            "Increases via VEGF-driven angiogenesis."
        ),
        node_type=NodeType.MEDIATOR,
        units="vessels/mm²",
        observable=True,
    ))

    dag.add_node(CausalNode(
        name="Survival",
        description=(
            "Cell survival probability (inverse of apoptosis rate). "
            "Enhanced by EGFR signaling through PI3K/AKT pathway."
        ),
        node_type=NodeType.OUTCOME,
        units="probability",
        observable=False,  # Inferred from cell death events
    ))

    # =========================================================================
    # EDGES (Causal Relationships)
    # =========================================================================

    # --- Step 1: ecDNA → EGFR expression (gene dosage) ---
    dag.add_edge(CausalEdge(
        source="ecDNA_EGFR",
        target="EGFR_expression",
        effect_type=EffectType.LINEAR,
        parameter_symbol="κ",
        parameter_name="EGFR_per_ecDNA_copy",
        parameter_value=egfr_per_copy,
        equation_template="EGFR = {base} + κ·ecDNA + noise",
        description=(
            "Gene dosage effect: each ecDNA copy adds κ units of EGFR expression. "
            "Includes transcriptional noise (log-normal). "
            "Ref: Hung et al. 2021 Nature (ecDNA forms transcriptional hubs)"
        ),
        sign="+",
    ))

    # --- Step 2: EGFR expression → Proliferation ---
    dag.add_edge(CausalEdge(
        source="EGFR_expression",
        target="Proliferation",
        effect_type=EffectType.LOG,
        parameter_symbol="α",
        parameter_name="ecDNA_effect_on_division",
        parameter_value=alpha,
        equation_template="T_div = T_base / (1 + α·log₂(1 + EGFR))",
        description=(
            "EGFR expression drives proliferation via RAS/MAPK "
            "and PI3K/AKT pathways. Saturating effect (log transform). "
            "Ref: Li et al. 2022 Nature Cell Biology"
        ),
        sign="+",
    ))

    # --- Proliferation depletes oxygen ---
    dag.add_edge(CausalEdge(
        source="Proliferation",
        target="O2_local",
        effect_type=EffectType.LINEAR,
        parameter_symbol="ω",
        parameter_name="O2_consumption_amol_hr",
        parameter_value=o2_consumption,
        equation_template="dO2/dt -= ω·Proliferation (consumption)",
        description=(
            "Proliferating cells consume oxygen for metabolism. Higher "
            "proliferation rate depletes local O2 faster."
        ),
        sign="-",
    ))

    # --- Low O2 creates hypoxia ---
    dag.add_edge(CausalEdge(
        source="O2_local",
        target="Hypoxia",
        effect_type=EffectType.THRESHOLD,
        parameter_symbol="τ",
        parameter_name="hypoxia_threshold_mmHg",
        parameter_value=hypoxia_threshold,
        equation_template="Hypoxia = 1 if O2 < τ else 0",
        description=(
            "Hypoxia defined as O2 below threshold. Triggers HIF-1α "
            "stabilization and downstream transcriptional program."
        ),
        sign="-",
    ))

    # --- EGFR expression drives invasion capacity ---
    dag.add_edge(CausalEdge(
        source="EGFR_expression",
        target="Invasion",
        effect_type=EffectType.MULTIPLICATIVE,
        parameter_symbol="δ",
        parameter_name="ecDNA_effect_on_migration",
        parameter_value=delta,
        equation_template="v_base_eff = v_base · (1 + δ·EGFR)",
        description=(
            "EGFR signaling enhances invasion capacity through EMT "
            "and MMP expression. This is the baseline invasion capacity "
            "before hypoxia modulation."
        ),
        sign="+",
    ))

    # --- Hypoxia drives invasion ("Go" response) ---
    dag.add_edge(CausalEdge(
        source="Hypoxia",
        target="Invasion",
        effect_type=EffectType.MULTIPLICATIVE,
        parameter_symbol="ψ",
        parameter_name="hypoxia_invasion_boost",
        parameter_value=2.0,
        equation_template="v_mig = v_base_eff · ψ (if hypoxic)",
        description=(
            "'Go or Grow' hypothesis: Hypoxia triggers invasion as escape "
            "response. HIF-1α activates EMT, MMPs, and migration machinery. "
            "Multiplicative with EGFR effect. "
            "Ref: Kathagen-Buhmann et al. 2016 Neuro-Oncology"
        ),
        sign="+",
    ))

    # --- Hypoxia suppresses proliferation ("Grow" inhibited) ---
    dag.add_edge(CausalEdge(
        source="Hypoxia",
        target="Proliferation",
        effect_type=EffectType.INVERSE,
        parameter_symbol="η",
        parameter_name="hypoxia_prolif_suppression",
        parameter_value=0.5,
        equation_template="Proliferation *= 1/(1 + η·Hypoxia)",
        description=(
            "'Go or Grow' hypothesis: Hypoxia suppresses proliferation due "
            "to metabolic constraints and phenotype switch to invasion."
        ),
        sign="-",
    ))

    # --- EGFR expression enhances VEGF secretion capacity ---
    dag.add_edge(CausalEdge(
        source="EGFR_expression",
        target="VEGF",
        effect_type=EffectType.SQRT,
        parameter_symbol="β",
        parameter_name="ecDNA_effect_on_VEGF",
        parameter_value=beta,
        equation_template="VEGF_capacity = base · (1 + β·√EGFR)",
        description=(
            "EGFR signaling can enhance HIF-1α stability and VEGF transcription "
            "even under normoxia. This sets the maximum VEGF secretion capacity."
        ),
        sign="+",
    ))

    # --- Hypoxia induces full VEGF secretion ---
    dag.add_edge(CausalEdge(
        source="Hypoxia",
        target="VEGF",
        effect_type=EffectType.MULTIPLICATIVE,
        parameter_symbol="φ",
        parameter_name="hypoxia_VEGF_multiplier",
        parameter_value=5.0,  # Full rate vs 20% basal
        equation_template="VEGF_secr = VEGF_capacity · (0.2 + 0.8·Hypoxia)",
        description=(
            "HIF-1α is stabilized under hypoxia and fully activates VEGF "
            "transcription. Under normoxia, only 20% basal secretion. "
            "Hypoxia multiplies by ~5x (full rate)."
        ),
        sign="+",
    ))

    # --- VEGF drives angiogenesis ---
    dag.add_edge(CausalEdge(
        source="VEGF",
        target="Vasculature",
        effect_type=EffectType.THRESHOLD,
        parameter_symbol="θ",
        parameter_name="angiogenesis_threshold_nM",
        parameter_value=angio_threshold,
        equation_template="Sprouting if VEGF > θ",
        description=(
            "VEGF above threshold triggers endothelial sprouting and "
            "new vessel formation (angiogenesis)."
        ),
        sign="+",
    ))

    # --- Vasculature delivers oxygen ---
    dag.add_edge(CausalEdge(
        source="Vasculature",
        target="O2_local",
        effect_type=EffectType.LINEAR,
        parameter_symbol="q",
        parameter_name="q_O2_transfer_per_hr",
        parameter_value=q_o2,
        equation_template="dO2/dt += q·Vessels·(O2_blood - O2_local)",
        description=(
            "Blood vessels deliver oxygen. More vessels = better "
            "oxygenation. Creates negative feedback with proliferation."
        ),
        sign="+",
    ))

    # --- EGFR expression enhances survival (PI3K/AKT) ---
    dag.add_edge(CausalEdge(
        source="EGFR_expression",
        target="Survival",
        effect_type=EffectType.LOG,
        parameter_symbol="γ",
        parameter_name="ecDNA_effect_on_survival",
        parameter_value=gamma,
        equation_template="Apoptosis_rate = base / (1 + γ·log₂(1 + EGFR))",
        description=(
            "EGFR signaling activates PI3K/AKT pathway which phosphorylates and "
            "inhibits pro-apoptotic proteins (BAD, FOXO)."
        ),
        sign="+",
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
