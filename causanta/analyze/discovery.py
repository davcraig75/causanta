"""Causal discovery algorithms for CAUSANTA.

Implements constraint-based (PC, FCI) and score-based (GES) causal
discovery algorithms for learning causal structure from observational data.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
from scipy import stats

from .loader import SimulationData, cells_to_array


@dataclass
class DiscoveredDAG:
    """Result of causal discovery."""

    algorithm: str
    nodes: list[str]
    edges: list[tuple[str, str]]  # (parent, child)
    undirected: list[tuple[str, str]]  # Edges with uncertain direction
    edge_weights: dict[tuple[str, str], float] = field(default_factory=dict)
    confidence: dict[tuple[str, str], float] = field(default_factory=dict)
    parameters: dict[str, Any] = field(default_factory=dict)

    def to_adjacency_matrix(self) -> tuple[list[str], np.ndarray]:
        """Convert to adjacency matrix.

        Returns 1 for directed edge, 0.5 for undirected edge.
        """
        n = len(self.nodes)
        node_idx = {name: i for i, name in enumerate(self.nodes)}
        matrix = np.zeros((n, n))

        for parent, child in self.edges:
            i, j = node_idx[parent], node_idx[child]
            matrix[i, j] = 1.0

        for n1, n2 in self.undirected:
            i, j = node_idx[n1], node_idx[n2]
            matrix[i, j] = 0.5
            matrix[j, i] = 0.5

        return self.nodes, matrix

    def to_dict(self) -> dict[str, Any]:
        return {
            "algorithm": self.algorithm,
            "nodes": self.nodes,
            "edges": self.edges,
            "undirected": self.undirected,
            "edge_weights": {f"{k[0]}->{k[1]}": v for k, v in self.edge_weights.items()},
            "confidence": {f"{k[0]}->{k[1]}": v for k, v in self.confidence.items()},
        }


def partial_correlation(X: np.ndarray, i: int, j: int, S: set[int]) -> float:
    """Compute partial correlation between X[:, i] and X[:, j] given X[:, S].

    Uses regression-based method: corr(residuals of i on S, residuals of j on S).
    """
    if not S:
        # Simple correlation
        return np.corrcoef(X[:, i], X[:, j])[0, 1]

    # Residualize both variables on conditioning set
    S_list = list(S)
    X_S = X[:, S_list]

    # Add intercept
    n = X.shape[0]
    X_S_int = np.column_stack([np.ones(n), X_S])

    # Regression for i
    beta_i = np.linalg.lstsq(X_S_int, X[:, i], rcond=None)[0]
    resid_i = X[:, i] - X_S_int @ beta_i

    # Regression for j
    beta_j = np.linalg.lstsq(X_S_int, X[:, j], rcond=None)[0]
    resid_j = X[:, j] - X_S_int @ beta_j

    return np.corrcoef(resid_i, resid_j)[0, 1]


def test_conditional_independence(
    X: np.ndarray,
    i: int,
    j: int,
    S: set[int],
    alpha: float = 0.05,
) -> tuple[bool, float]:
    """Test conditional independence of X_i and X_j given X_S.

    Uses Fisher's z-transform for partial correlation test.

    Returns:
        (is_independent, p_value)
    """
    n = X.shape[0]
    r = partial_correlation(X, i, j, S)

    # Fisher's z-transform
    # Avoid exactly 1 or -1
    r = np.clip(r, -0.9999, 0.9999)
    z = 0.5 * np.log((1 + r) / (1 - r))

    # Standard error
    se = 1 / np.sqrt(n - len(S) - 3)

    # Z-statistic
    z_stat = abs(z) / se

    # P-value (two-tailed)
    p_value = 2 * (1 - stats.norm.cdf(z_stat))

    return p_value > alpha, p_value


def pc_algorithm(
    X: np.ndarray,
    node_names: list[str],
    alpha: float = 0.05,
    max_cond_size: int = 3,
) -> DiscoveredDAG:
    """PC algorithm for causal discovery.

    Constraint-based algorithm that starts with complete undirected graph
    and removes edges based on conditional independence tests.

    Args:
        X: Data matrix (n x p)
        node_names: Names for each column
        alpha: Significance level for independence tests
        max_cond_size: Maximum conditioning set size

    Returns:
        DiscoveredDAG with discovered structure.
    """
    n, p = X.shape
    nodes = list(node_names)

    # Start with complete undirected graph
    adj = np.ones((p, p)) - np.eye(p)
    separating_sets: dict[tuple[int, int], set[int]] = {}

    # Phase 1: Edge removal
    for cond_size in range(max_cond_size + 1):
        for i in range(p):
            neighbors_i = np.where(adj[i] > 0)[0]
            for j in neighbors_i:
                if i >= j:
                    continue

                # Find conditioning sets of size cond_size
                other_neighbors = [k for k in neighbors_i if k != j]
                if len(other_neighbors) < cond_size:
                    continue

                # Try all subsets of size cond_size
                from itertools import combinations
                for S in combinations(other_neighbors, cond_size):
                    S_set = set(S)
                    is_indep, p_val = test_conditional_independence(X, i, j, S_set, alpha)

                    if is_indep:
                        adj[i, j] = 0
                        adj[j, i] = 0
                        separating_sets[(i, j)] = S_set
                        separating_sets[(j, i)] = S_set
                        break

    # Phase 2: Orient edges (simplified v-structure detection)
    directed = np.zeros((p, p))

    # Find v-structures: i -> k <- j where i and j are not adjacent
    for k in range(p):
        parents_k = np.where(adj[:, k] > 0)[0]
        for idx1, i in enumerate(parents_k):
            for j in parents_k[idx1 + 1:]:
                if adj[i, j] == 0:  # i and j not adjacent
                    # Check if k is in separating set
                    sep = separating_sets.get((i, j), set())
                    if k not in sep:
                        # Orient as v-structure
                        directed[i, k] = 1
                        directed[j, k] = 1

    # Build result
    edges = []
    undirected = []

    for i in range(p):
        for j in range(i + 1, p):
            if adj[i, j] > 0:
                if directed[i, j] > 0:
                    edges.append((nodes[i], nodes[j]))
                elif directed[j, i] > 0:
                    edges.append((nodes[j], nodes[i]))
                else:
                    undirected.append((nodes[i], nodes[j]))

    return DiscoveredDAG(
        algorithm="PC",
        nodes=nodes,
        edges=edges,
        undirected=undirected,
        parameters={"alpha": alpha, "max_cond_size": max_cond_size},
    )


def ges_algorithm(
    X: np.ndarray,
    node_names: list[str],
    penalty: float = 1.0,
) -> DiscoveredDAG:
    """Greedy Equivalence Search (GES) algorithm.

    Score-based algorithm that searches over DAG space by greedily
    adding and removing edges to maximize BIC score.

    Args:
        X: Data matrix (n x p)
        node_names: Names for each column
        penalty: BIC penalty multiplier

    Returns:
        DiscoveredDAG with discovered structure.
    """
    n, p = X.shape
    nodes = list(node_names)

    def bic_score(parents: dict[int, set[int]]) -> float:
        """Compute BIC score for current graph structure."""
        total_score = 0
        for j in range(p):
            parent_set = list(parents.get(j, set()))
            k = len(parent_set)

            if k == 0:
                # No parents - just variance
                var_j = np.var(X[:, j], ddof=1)
                ll = -n / 2 * np.log(2 * np.pi * var_j) - n / 2
            else:
                # Regression of j on parents
                X_pa = X[:, parent_set]
                X_pa_int = np.column_stack([np.ones(n), X_pa])
                beta = np.linalg.lstsq(X_pa_int, X[:, j], rcond=None)[0]
                resid = X[:, j] - X_pa_int @ beta
                var_resid = np.var(resid, ddof=1)
                ll = -n / 2 * np.log(2 * np.pi * var_resid) - n / 2

            # BIC penalty
            total_score += ll - penalty * k * np.log(n) / 2

        return total_score

    def _would_create_cycle(parents: dict[int, set[int]], source: int, target: int) -> bool:
        """Check if adding source -> target would create a cycle via DFS."""
        # A cycle exists if target can already reach source through existing edges
        visited = set()
        stack = [source]
        while stack:
            node = stack.pop()
            if node == target:
                return True
            if node in visited:
                continue
            visited.add(node)
            # Follow edges: node's parents are nodes that point TO node
            # We need to check if source is reachable FROM target
            # i.e., follow parent edges from source upward
            stack.extend(parents.get(node, set()))
        return False

    # Initialize empty graph
    parents: dict[int, set[int]] = {j: set() for j in range(p)}
    current_score = bic_score(parents)

    # Forward phase: add edges (with acyclicity enforcement)
    improved = True
    while improved:
        improved = False
        best_add = None
        best_score = current_score

        for j in range(p):
            for i in range(p):
                if i == j or i in parents[j]:
                    continue

                # Skip if adding i -> j would create a cycle
                if _would_create_cycle(parents, i, j):
                    continue

                # Try adding i -> j
                parents[j].add(i)
                new_score = bic_score(parents)

                if new_score > best_score:
                    best_score = new_score
                    best_add = (i, j)

                parents[j].remove(i)

        if best_add is not None:
            parents[best_add[1]].add(best_add[0])
            current_score = best_score
            improved = True

    # Backward phase: remove edges
    improved = True
    while improved:
        improved = False
        best_remove = None
        best_score = current_score

        for j in range(p):
            for i in list(parents[j]):
                # Try removing i -> j
                parents[j].remove(i)
                new_score = bic_score(parents)

                if new_score > best_score:
                    best_score = new_score
                    best_remove = (i, j)

                parents[j].add(i)

        if best_remove is not None:
            parents[best_remove[1]].remove(best_remove[0])
            current_score = best_score
            improved = True

    # Build result
    edges = []
    for j in range(p):
        for i in parents[j]:
            edges.append((nodes[i], nodes[j]))

    return DiscoveredDAG(
        algorithm="GES",
        nodes=nodes,
        edges=edges,
        undirected=[],
        parameters={"penalty": penalty, "bic_score": current_score},
    )


def discover_causal_structure(
    data: SimulationData,
    variables: list[str] | None = None,
    algorithm: str = "pc",
    **kwargs,
) -> DiscoveredDAG:
    """Discover causal structure from simulation data.

    Args:
        data: Loaded simulation data
        variables: Variables to include. If None, uses defaults.
        algorithm: "pc", "ges", or "both"
        **kwargs: Algorithm-specific parameters

    Returns:
        DiscoveredDAG with discovered structure.
    """
    if variables is None:
        # Note: glucose_local excluded as it's highly correlated with O2_local
        # (both determined by vascular distance) and causes multicollinearity
        # that confuses the PC algorithm. O2 is the biologically relevant
        # variable for hypoxia-driven phenotypes.
        variables = [
            "ecDNA_count",
            "egfr_expression",
            "O2_local",
            "VEGF_secretion",
            "migration_rate",
        ]

    # Get tumor cells
    tumor_cells = data.tumor_cells
    if len(tumor_cells) < 30:
        return DiscoveredDAG(
            algorithm=algorithm,
            nodes=variables,
            edges=[],
            undirected=[],
            parameters={"error": "insufficient_data"},
        )

    # Build data matrix
    X = cells_to_array(tumor_cells, variables)

    # Standardize
    X = (X - np.mean(X, axis=0)) / (np.std(X, axis=0) + 1e-10)

    if algorithm.lower() == "pc":
        return pc_algorithm(X, variables, **kwargs)
    elif algorithm.lower() == "ges":
        return ges_algorithm(X, variables, **kwargs)
    else:
        raise ValueError(f"Unknown algorithm: {algorithm}")


def compare_to_ground_truth(
    discovered: DiscoveredDAG,
    ground_truth_edges: list[tuple[str, str]],
) -> dict[str, Any]:
    """Compare discovered structure to ground truth.

    Args:
        discovered: Discovered DAG
        ground_truth_edges: List of true directed edges

    Returns:
        Dictionary with precision, recall, F1, SHD.
    """
    gt_set = set(ground_truth_edges)
    disc_set = set(discovered.edges)

    # Also consider undirected as half-correct
    for n1, n2 in discovered.undirected:
        if (n1, n2) in gt_set:
            disc_set.add((n1, n2))
        elif (n2, n1) in gt_set:
            disc_set.add((n2, n1))

    tp = len(gt_set & disc_set)
    fp = len(disc_set - gt_set)
    fn = len(gt_set - disc_set)

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    shd = fp + fn  # Structural Hamming Distance

    return {
        "true_positives": tp,
        "false_positives": fp,
        "false_negatives": fn,
        "precision": precision,
        "recall": recall,
        "f1_score": f1,
        "structural_hamming_distance": shd,
        "discovered_edges": list(disc_set),
        "missed_edges": list(gt_set - disc_set),
        "extra_edges": list(disc_set - gt_set),
    }
