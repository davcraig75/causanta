"""Environment fields and diffusion solver for CAUSANTA.

Implements reaction-diffusion PDEs on the environment grid using
implicit Locally One-Dimensional (LOD) operator splitting with
the vectorized Thomas algorithm for unconditional stability.

Performance: Uses Numba JIT compilation when available for ~10x speedup.
"""

from __future__ import annotations

import numpy as np

from .config import EnvironmentConfig, SimulationConfig
from .domain import Domain

# Import Numba-accelerated kernels if available
try:
    from ._numba_kernels import (
        HAS_NUMBA,
        get_thomas_solver,
        thomas_solve_batch_numba,
        thomas_solve_batch_parallel,
    )
except ImportError:
    HAS_NUMBA = False
    get_thomas_solver = None


class EnvironmentFields:
    """Stores all environment field arrays on the grid."""

    def __init__(self, domain: Domain, config: EnvironmentConfig):
        self.domain = domain
        ny, nx = domain.ny, domain.nx

        # Primary diffusing fields (shape: ny, nx) — field[row, col]
        sub = config.substrates
        self.O2 = np.full((ny, nx), sub["O2"].initial_value, dtype=np.float64)
        self.glucose = np.full((ny, nx), sub["glucose"].initial_value, dtype=np.float64)
        self.VEGF = np.full((ny, nx), sub["VEGF"].initial_value, dtype=np.float64)
        self.lactate = np.full((ny, nx), sub["lactate"].initial_value, dtype=np.float64)

        # Static fields
        self.ECM_density = np.full((ny, nx), 0.5, dtype=np.float64)
        self.vascular_density = np.zeros((ny, nx), dtype=np.float64)

        # Derived
        self.pH = np.full((ny, nx), 7.4, dtype=np.float64)

        # Occupancy tracking
        self.is_occupied = np.zeros((ny, nx), dtype=bool)
        self.occupant_cell_id = np.full((ny, nx), -1, dtype=np.int64)
        self.occupant_type = np.full((ny, nx), -1, dtype=np.int8)

    def get_field(self, name: str) -> np.ndarray:
        """Get a field array by name."""
        return getattr(self, name)

    def update_derived(self) -> None:
        """Recompute derived fields from primary fields."""
        # pH = 7.4 - 0.02 * lactate
        np.subtract(7.4, 0.02 * self.lactate, out=self.pH)

    def clear_occupancy(self) -> None:
        """Reset occupancy arrays."""
        self.is_occupied[:] = False
        self.occupant_cell_id[:] = -1
        self.occupant_type[:] = -1


def _thomas_solve_batch_numpy(
    a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray
) -> np.ndarray:
    """Pure NumPy Thomas algorithm for batch tridiagonal systems.

    Fallback when Numba is not available.

    Args:
        a: Sub-diagonal coefficients, shape (M, N-1)
        b: Main diagonal coefficients, shape (M, N)
        c: Super-diagonal coefficients, shape (M, N-1)
        d: Right-hand side, shape (M, N)

    Returns:
        Solution x, shape (M, N)
    """
    n = d.shape[1]

    # Work on copies to avoid modifying inputs
    c_p = np.empty_like(c)
    d_p = np.empty_like(d)

    # Forward sweep
    c_p[:, 0] = c[:, 0] / b[:, 0]
    d_p[:, 0] = d[:, 0] / b[:, 0]

    for i in range(1, n):
        denom = b[:, i] - a[:, i - 1] * c_p[:, i - 1]
        if i < n - 1:
            c_p[:, i] = c[:, i] / denom
        d_p[:, i] = (d[:, i] - a[:, i - 1] * d_p[:, i - 1]) / denom

    # Back substitution
    x = np.empty_like(d)
    x[:, -1] = d_p[:, -1]
    for i in range(n - 2, -1, -1):
        x[:, i] = d_p[:, i] - c_p[:, i] * x[:, i + 1]

    return x


def _thomas_solve_batch(
    a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray
) -> np.ndarray:
    """Dispatch to fastest available Thomas solver.

    Uses Numba JIT-compiled version when available (~10x faster),
    falls back to pure NumPy implementation otherwise.
    """
    if HAS_NUMBA and get_thomas_solver is not None:
        solver = get_thomas_solver(d.shape[0])
        return solver(a, b, c, d)
    return _thomas_solve_batch_numpy(a, b, c, d)


class DiffusionSolver:
    """Implicit LOD diffusion solver with vectorized Thomas algorithm."""

    def __init__(self, domain: Domain, config: EnvironmentConfig):
        self.domain = domain
        self.dx = float(domain.env_grid_um)
        self.config = config
        self._substrate_names = ["O2", "glucose", "VEGF", "lactate"]

    def _sweep_x(self, field: np.ndarray, r: float, boundary_value: float) -> np.ndarray:
        """Implicit sweep in x-direction for all rows simultaneously.

        Solves: -r * u_{i-1}^{n+1} + (1+2r) * u_i^{n+1} - r * u_{i+1}^{n+1} = u_i^n
        with Dirichlet BCs at i=0 and i=nx-1.
        """
        ny, nx = field.shape
        if nx <= 2:
            return field.copy()

        # Interior points only (indices 1..nx-2)
        n_interior = nx - 2
        n_rows = ny

        # Build tridiagonal coefficients for interior points
        # Main diagonal: (1 + 2r)
        b = np.full((n_rows, n_interior), 1.0 + 2.0 * r)
        # Sub/super diagonal: -r
        a = np.full((n_rows, n_interior - 1), -r)
        c = np.full((n_rows, n_interior - 1), -r)

        # RHS = current field values at interior points
        d = field[:, 1:-1].copy()
        # Add BC contributions: left boundary affects first interior point
        d[:, 0] += r * boundary_value
        # Right boundary affects last interior point
        d[:, -1] += r * boundary_value

        # Solve
        result = field.copy()
        result[:, 0] = boundary_value
        result[:, -1] = boundary_value
        if n_interior > 0:
            result[:, 1:-1] = _thomas_solve_batch(a, b, c, d)
        return result

    def _sweep_y(self, field: np.ndarray, r: float, boundary_value: float) -> np.ndarray:
        """Implicit sweep in y-direction for all columns simultaneously.

        Same as sweep_x but transposed.
        """
        ny, nx = field.shape
        if ny <= 2:
            return field.copy()

        n_interior = ny - 2
        n_cols = nx

        # Transpose to work column-wise: shape (nx, ny) -> solve along dim 1
        ft = field.T  # shape (nx, ny)

        b = np.full((n_cols, n_interior), 1.0 + 2.0 * r)
        a = np.full((n_cols, n_interior - 1), -r)
        c = np.full((n_cols, n_interior - 1), -r)

        d = ft[:, 1:-1].copy()
        d[:, 0] += r * boundary_value
        d[:, -1] += r * boundary_value

        result_t = ft.copy()
        result_t[:, 0] = boundary_value
        result_t[:, -1] = boundary_value
        if n_interior > 0:
            result_t[:, 1:-1] = _thomas_solve_batch(a, b, c, d)
        return result_t.T

    def solve_substep(
        self,
        field: np.ndarray,
        substrate_name: str,
        dt: float,
        source: np.ndarray,
        sink: np.ndarray,
    ) -> np.ndarray:
        """One LOD diffusion sub-step for a single substrate.

        1. X-sweep (implicit, half the diffusion)
        2. Y-sweep (implicit, half the diffusion)
        3. Apply source/sink and decay (operator split)
        4. Enforce BCs and clamp
        """
        sub = self.config.substrates[substrate_name]
        D = sub.diffusion_coeff_um2_hr
        lam = sub.decay_rate_per_hr
        bv = sub.boundary_value

        r = D * dt / (self.dx * self.dx)

        # LOD: split diffusion into x and y half-steps
        field = self._sweep_x(field, r / 2.0, bv)
        field = self._sweep_y(field, r / 2.0, bv)

        # Reaction terms (explicit operator split)
        field += dt * (source - sink - lam * field)

        # Enforce boundary conditions
        field[0, :] = bv
        field[-1, :] = bv
        field[:, 0] = bv
        field[:, -1] = bv

        # Clamp to physical range
        np.clip(field, sub.clamp_min, sub.clamp_max, out=field)

        return field

    def step(
        self,
        env: EnvironmentFields,
        sources: dict[str, np.ndarray],
        sinks: dict[str, np.ndarray],
        dt_diffusion: float,
        n_substeps: int,
    ) -> None:
        """Advance all substrate fields by one main time step (1 hour).

        Args:
            env: Environment fields to update in-place
            sources: Source arrays per substrate
            sinks: Sink arrays per substrate
            dt_diffusion: Diffusion sub-step size (hours)
            n_substeps: Number of sub-steps per hour
        """
        for _ in range(n_substeps):
            for name in self._substrate_names:
                field = env.get_field(name)
                src = sources.get(name, np.zeros_like(field))
                snk = sinks.get(name, np.zeros_like(field))
                updated = self.solve_substep(field, name, dt_diffusion, src, snk)
                # Write back in-place
                field[:] = updated

        env.update_derived()


def accumulate_sources_and_sinks(
    population,  # CellPopulation — import avoided for circular dependency
    env: EnvironmentFields,
    domain: Domain,
    config: SimulationConfig,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    """Compute source and sink arrays from cell population and vasculature.

    Called once per main time step; arrays are reused across all diffusion sub-steps.
    """
    ny, nx = domain.ny, domain.nx
    voxel_area = domain.env_grid_um ** 2

    # Initialize arrays
    O2_source = np.zeros((ny, nx))
    O2_sink = np.zeros((ny, nx))
    glucose_source = np.zeros((ny, nx))
    glucose_sink = np.zeros((ny, nx))
    VEGF_source = np.zeros((ny, nx))
    lactate_source = np.zeros((ny, nx))

    # NOTE ON UNITS: The reaction terms below mix physical units (mmHg, mM,
    # amol/hr) without explicit conversion factors. The source/sink terms are
    # phenomenological — parameter values are tuned for qualitative behavior,
    # not dimensional consistency. This is acceptable for a causal benchmark
    # but means parameters cannot be interpreted as literal physical quantities.

    # Vascular supply (vectorized over grid)
    O2_source += config.environment.q_O2_transfer_per_hr * env.vascular_density * (
        config.environment.O2_blood_mmHg - env.O2
    )
    glucose_source += config.environment.q_glucose_transfer_per_hr * env.vascular_density * (
        5.0 - env.glucose  # glucose blood level
    )

    # Cell contributions
    hypoxia_thresh = config.environment.hypoxia_threshold_mmHg
    for cell in population.iter_living():
        gi, gj = domain.um_to_grid(cell.x, cell.y)
        tp = config.cell_types.get(cell.cell_type)
        if tp is None:
            continue

        # Oxygen consumption
        O2_sink[gj, gi] += tp.O2_consumption_amol_hr / voxel_area
        # Glucose consumption
        glucose_sink[gj, gi] += tp.glucose_consumption_amol_hr / voxel_area

        # VEGF secretion: cell.VEGF_secretion now already includes the
        # hypoxia gating (0.2x normoxic, 1.0x hypoxic) from
        # modulate_vegf_secretion in ecdna.py. Just pass through.
        if cell.VEGF_secretion > 0:
            VEGF_source[gj, gi] += cell.VEGF_secretion / voxel_area

        # Lactate production
        if tp.lactate_production_amol_hr > 0:
            lactate_source[gj, gi] += tp.lactate_production_amol_hr / voxel_area

    sources = {
        "O2": O2_source,
        "glucose": glucose_source,
        "VEGF": VEGF_source,
        "lactate": lactate_source,
    }
    sinks = {
        "O2": O2_sink,
        "glucose": glucose_sink,
        "VEGF": np.zeros((ny, nx)),
        "lactate": np.zeros((ny, nx)),
    }
    return sources, sinks
