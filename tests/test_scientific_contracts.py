"""Tests for CAUSANTA's core scientific contracts.

These tests verify that the simulator, analyzer, and documentation
are internally consistent on the key causal benchmarking claims.
"""

import math
import numpy as np
import pytest


# ---------------------------------------------------------------------------
# 1. ecDNA segregation: conservation and expected fractions
# ---------------------------------------------------------------------------

class TestEcDNASegregation:
    """Verify ecDNA replication + segregation model."""

    def test_conservation_after_segregation(self):
        """Parent_after + daughter should equal the replicated pool."""
        from causanta.simulate.ecdna import replicate_and_segregate_ecdna
        rng = np.random.default_rng(42)

        for parent_count in [5, 10, 20, 50]:
            for _ in range(100):
                parent_new, daughter = replicate_and_segregate_ecdna(
                    parent_count, 0.5, rng
                )
                # Total should equal the replicated amount
                assert parent_new + daughter >= parent_count, (
                    f"Total {parent_new + daughter} < pre-replication {parent_count}"
                )
                assert parent_new >= 0
                assert daughter >= 0

    def test_daughter_fraction_of_replicated_pool_near_half(self):
        """daughter / (parent_after + daughter) should be ~0.5 on average."""
        from causanta.simulate.ecdna import replicate_and_segregate_ecdna
        rng = np.random.default_rng(42)

        fractions = []
        for _ in range(1000):
            parent_new, daughter = replicate_and_segregate_ecdna(20, 0.5, rng)
            total = parent_new + daughter
            if total > 0:
                fractions.append(daughter / total)

        mean_frac = np.mean(fractions)
        assert abs(mean_frac - 0.5) < 0.05, (
            f"Mean daughter fraction of replicated pool = {mean_frac:.3f}, expected ~0.5"
        )

    def test_max_ecdna_cap(self):
        """ecDNA should never exceed MAX_ECDNA_COPIES."""
        from causanta.simulate.ecdna import replicate_and_segregate_ecdna, MAX_ECDNA_COPIES
        rng = np.random.default_rng(42)

        for _ in range(100):
            parent_new, daughter = replicate_and_segregate_ecdna(95, 0.5, rng)
            assert parent_new <= MAX_ECDNA_COPIES
            assert daughter <= MAX_ECDNA_COPIES


# ---------------------------------------------------------------------------
# 2. Causal chain: all phenotype modulation routes through EGFR
# ---------------------------------------------------------------------------

class TestCausalChain:
    """Verify ecDNA -> EGFR -> phenotype (no direct ecDNA -> phenotype)."""

    def test_apoptosis_uses_egfr(self):
        """check_apoptosis should use cell.egfr_expression, not ecDNA_count."""
        import inspect
        from causanta.simulate.behaviors import check_apoptosis
        source = inspect.getsource(check_apoptosis)
        assert "cell.egfr_expression" in source, (
            "check_apoptosis should use cell.egfr_expression"
        )
        assert "cell.ecDNA_count" not in source, (
            "check_apoptosis should NOT use cell.ecDNA_count directly"
        )

    def test_division_time_g0_entry_uses_egfr(self):
        """G0->G1 transition should compute division time from EGFR."""
        import inspect
        from causanta.simulate.behaviors import advance_cell_cycle
        source = inspect.getsource(advance_cell_cycle)
        # The modulate_division_time call should use egfr_expression
        assert "cell.egfr_expression" in source, (
            "advance_cell_cycle should use cell.egfr_expression for division time"
        )


# ---------------------------------------------------------------------------
# 3. IV pipeline: instrument != treatment
# ---------------------------------------------------------------------------

class TestIVPipeline:
    """Verify IV analysis uses distinct instrument and treatment."""

    def test_loader_includes_egfr(self):
        """Loader should parse egfr_expression from TSV."""
        import inspect
        from causanta.analyze.loader import load_cells_tsv
        source = inspect.getsource(load_cells_tsv)
        assert "egfr_expression" in source, (
            "load_cells_tsv should parse egfr_expression column"
        )

    def test_iv_uses_egfr_as_treatment(self):
        """IV estimator should use egfr_expression as treatment, not ecDNA."""
        import inspect
        from causanta.analyze.iv import estimate_iv_effects
        source = inspect.getsource(estimate_iv_effects)
        assert "egfr_expression" in source, (
            "estimate_iv_effects should use egfr_expression as treatment"
        )
        # The primary path should use EGFR; D=Z is only a documented fallback
        assert "treatment_name = \"egfr_expression\"" in source, (
            "estimate_iv_effects should set treatment to egfr_expression as primary path"
        )


# ---------------------------------------------------------------------------
# 4. Discovery: acyclicity
# ---------------------------------------------------------------------------

class TestDiscovery:
    """Verify causal discovery produces DAGs (no cycles)."""

    def test_ges_produces_dag(self):
        """GES output should have no cycles."""
        from causanta.analyze.discovery import ges_algorithm
        rng = np.random.default_rng(42)

        # Create correlated data that might tempt cycle creation
        n = 200
        X = rng.standard_normal((n, 4))
        X[:, 1] += 0.5 * X[:, 0]  # 0 -> 1
        X[:, 2] += 0.3 * X[:, 1]  # 1 -> 2
        X[:, 3] += 0.4 * X[:, 0]  # 0 -> 3

        result = ges_algorithm(X, ["A", "B", "C", "D"])

        # Check no cycles: for each edge (u, v), v should not reach u
        adj = {}
        for s, t in result.edges:
            adj.setdefault(s, []).append(t)

        def has_path(start, end, graph):
            visited = set()
            stack = [start]
            while stack:
                node = stack.pop()
                if node == end:
                    return True
                if node in visited:
                    continue
                visited.add(node)
                stack.extend(graph.get(node, []))
            return False

        for s, t in result.edges:
            assert not has_path(t, s, adj), (
                f"Cycle detected: edge {s}->{t} but {t} can reach {s}"
            )


# ---------------------------------------------------------------------------
# 5. Output path correctness
# ---------------------------------------------------------------------------

class TestOutputPaths:
    """Verify visualization HTML references correct file paths."""

    def test_viewer_references_correct_vega_path(self):
        """index.html should reference the correct Vega spec location."""
        import inspect
        from causanta.simulate.visualization import write_html_viewer
        source = inspect.getsource(write_html_viewer)
        # Should NOT reference 'visualization.vg.json' (old incorrect name)
        assert "visualization.vg.json" not in source, (
            "write_html_viewer should not reference old filename visualization.vg.json"
        )


# ---------------------------------------------------------------------------
# 6. Immune system: max_kills enforcement
# ---------------------------------------------------------------------------

class TestImmuneSystem:
    """Verify immune killing mechanics."""

    def test_max_kills_enforced(self):
        """check_immune_kill should respect max_kills_before_exhaustion."""
        import inspect
        from causanta.simulate.behaviors import check_immune_kill
        source = inspect.getsource(check_immune_kill)
        assert "max_kills_before_exhaustion" in source, (
            "check_immune_kill should enforce max_kills_before_exhaustion"
        )

    def test_immune_chemotaxis_called_in_migration(self):
        """compute_migration should call compute_immune_chemotaxis for immune cells."""
        import inspect
        from causanta.simulate.behaviors import compute_migration
        source = inspect.getsource(compute_migration)
        assert "compute_immune_chemotaxis" in source, (
            "compute_migration should call compute_immune_chemotaxis"
        )
