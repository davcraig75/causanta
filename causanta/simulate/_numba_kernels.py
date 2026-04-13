"""Numba-accelerated kernels for CAUSANTA simulation.

These JIT-compiled functions replace the pure Python implementations
for performance-critical inner loops.
"""

from __future__ import annotations

import numpy as np

# Try to import numba, fall back to pure Python if not available
try:
    from numba import njit, prange
    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False

    # Create no-op decorator that mimics njit signature
    def njit(*args, **kwargs):
        def decorator(func):
            return func
        if len(args) == 1 and callable(args[0]):
            return args[0]
        return decorator

    def prange(*args, **kwargs):
        return range(*args)


# ============================================================================
# Thomas Algorithm (Tridiagonal Solver)
# ============================================================================

@njit(cache=True)
def thomas_solve_batch_numba(
    a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray
) -> np.ndarray:
    """Numba-accelerated Thomas algorithm for batch tridiagonal systems.

    Solves M independent tridiagonal systems simultaneously.

    Args:
        a: Sub-diagonal coefficients, shape (M, N-1)
        b: Main diagonal coefficients, shape (M, N)
        c: Super-diagonal coefficients, shape (M, N-1)
        d: Right-hand side, shape (M, N)

    Returns:
        Solution x, shape (M, N)
    """
    m, n = d.shape

    # Work arrays
    c_prime = np.empty((m, n - 1), dtype=np.float64)
    d_prime = np.empty((m, n), dtype=np.float64)
    x = np.empty((m, n), dtype=np.float64)

    # Forward sweep - vectorized over rows
    for row in range(m):
        c_prime[row, 0] = c[row, 0] / b[row, 0]
        d_prime[row, 0] = d[row, 0] / b[row, 0]

        for i in range(1, n):
            denom = b[row, i] - a[row, i - 1] * c_prime[row, i - 1]
            if i < n - 1:
                c_prime[row, i] = c[row, i] / denom
            d_prime[row, i] = (d[row, i] - a[row, i - 1] * d_prime[row, i - 1]) / denom

        # Back substitution
        x[row, n - 1] = d_prime[row, n - 1]
        for i in range(n - 2, -1, -1):
            x[row, i] = d_prime[row, i] - c_prime[row, i] * x[row, i + 1]

    return x


@njit(parallel=True, cache=True)
def thomas_solve_batch_parallel(
    a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray
) -> np.ndarray:
    """Parallel Numba Thomas algorithm - each row solved independently.

    For large grids (>50 rows), parallel version is faster.
    """
    m, n = d.shape

    # Pre-allocate output
    x = np.empty((m, n), dtype=np.float64)

    # Each row is independent - parallelize
    for row in prange(m):
        # Local work arrays
        c_prime = np.empty(n - 1, dtype=np.float64)
        d_prime = np.empty(n, dtype=np.float64)

        # Forward sweep
        c_prime[0] = c[row, 0] / b[row, 0]
        d_prime[0] = d[row, 0] / b[row, 0]

        for i in range(1, n):
            denom = b[row, i] - a[row, i - 1] * c_prime[i - 1]
            if i < n - 1:
                c_prime[i] = c[row, i] / denom
            d_prime[i] = (d[row, i] - a[row, i - 1] * d_prime[i - 1]) / denom

        # Back substitution
        x[row, n - 1] = d_prime[n - 1]
        for i in range(n - 2, -1, -1):
            x[row, i] = d_prime[i] - c_prime[i] * x[row, i + 1]

    return x


# ============================================================================
# Diffusion Sweep Operations
# ============================================================================

@njit(cache=True)
def apply_reaction_terms(
    field: np.ndarray,
    source: np.ndarray,
    sink: np.ndarray,
    decay_rate: float,
    dt: float,
) -> np.ndarray:
    """Apply reaction terms (source, sink, decay) in-place.

    field += dt * (source - sink - decay_rate * field)
    """
    ny, nx = field.shape
    result = np.empty((ny, nx), dtype=np.float64)
    for j in range(ny):
        for i in range(nx):
            result[j, i] = field[j, i] + dt * (
                source[j, i] - sink[j, i] - decay_rate * field[j, i]
            )
    return result


@njit(cache=True)
def clamp_field(field: np.ndarray, min_val: float, max_val: float) -> None:
    """Clamp field values in-place."""
    ny, nx = field.shape
    for j in range(ny):
        for i in range(nx):
            if field[j, i] < min_val:
                field[j, i] = min_val
            elif field[j, i] > max_val:
                field[j, i] = max_val


# ============================================================================
# Cell Source/Sink Accumulation
# ============================================================================

@njit(cache=True)
def accumulate_cell_contributions(
    cell_x: np.ndarray,
    cell_y: np.ndarray,
    cell_types: np.ndarray,
    cell_is_alive: np.ndarray,
    cell_O2_local: np.ndarray,
    cell_VEGF_secretion: np.ndarray,
    type_O2_consumption: np.ndarray,
    type_glucose_consumption: np.ndarray,
    type_lactate_production: np.ndarray,
    env_grid_um: float,
    hypoxia_threshold: float,
    O2_sink: np.ndarray,
    glucose_sink: np.ndarray,
    VEGF_source: np.ndarray,
    lactate_source: np.ndarray,
) -> None:
    """Accumulate cell contributions to source/sink arrays.

    All arrays are modified in-place.
    """
    voxel_area = env_grid_um * env_grid_um
    n_cells = len(cell_x)
    ny, nx = O2_sink.shape

    for i in range(n_cells):
        if not cell_is_alive[i]:
            continue

        # Convert position to grid indices
        gi = int(cell_x[i] / env_grid_um)
        gj = int(cell_y[i] / env_grid_um)

        # Clamp to grid bounds
        if gi < 0:
            gi = 0
        elif gi >= nx:
            gi = nx - 1
        if gj < 0:
            gj = 0
        elif gj >= ny:
            gj = ny - 1

        cell_type = cell_types[i]

        # Accumulate contributions
        O2_sink[gj, gi] += type_O2_consumption[cell_type] / voxel_area
        glucose_sink[gj, gi] += type_glucose_consumption[cell_type] / voxel_area

        # VEGF secretion (hypoxia-dependent)
        if cell_VEGF_secretion[i] > 0 and cell_O2_local[i] < hypoxia_threshold:
            VEGF_source[gj, gi] += cell_VEGF_secretion[i] / voxel_area

        # Lactate production
        if type_lactate_production[cell_type] > 0:
            lactate_source[gj, gi] += type_lactate_production[cell_type] / voxel_area


# ============================================================================
# Gradient Computation
# ============================================================================

@njit(cache=True)
def compute_gradient_at_point(
    field: np.ndarray,
    x: float,
    y: float,
    grid_spacing: float,
) -> tuple:
    """Compute gradient of field at (x, y) using central differences.

    Returns (grad_x, grad_y).
    """
    ny, nx = field.shape

    # Convert to grid coordinates
    gx = x / grid_spacing
    gy = y / grid_spacing

    # Integer grid indices
    gi = int(gx)
    gj = int(gy)

    # Clamp to valid range for gradient computation
    gi = max(1, min(gi, nx - 2))
    gj = max(1, min(gj, ny - 2))

    # Central difference gradient
    grad_x = (field[gj, gi + 1] - field[gj, gi - 1]) / (2.0 * grid_spacing)
    grad_y = (field[gj + 1, gi] - field[gj - 1, gi]) / (2.0 * grid_spacing)

    return grad_x, grad_y


@njit(cache=True)
def bilinear_interpolate(
    field: np.ndarray,
    x: float,
    y: float,
    grid_spacing: float,
) -> float:
    """Bilinear interpolation of field at (x, y)."""
    ny, nx = field.shape

    # Convert to grid coordinates
    gx = x / grid_spacing
    gy = y / grid_spacing

    # Integer grid indices (lower-left corner)
    gi0 = int(gx)
    gj0 = int(gy)

    # Clamp to valid range
    gi0 = max(0, min(gi0, nx - 2))
    gj0 = max(0, min(gj0, ny - 2))

    gi1 = gi0 + 1
    gj1 = gj0 + 1

    # Fractional parts
    fx = gx - gi0
    fy = gy - gj0

    # Clamp fractions
    if fx < 0:
        fx = 0.0
    elif fx > 1:
        fx = 1.0
    if fy < 0:
        fy = 0.0
    elif fy > 1:
        fy = 1.0

    # Bilinear interpolation
    v00 = field[gj0, gi0]
    v10 = field[gj0, gi1]
    v01 = field[gj1, gi0]
    v11 = field[gj1, gi1]

    return (
        v00 * (1 - fx) * (1 - fy) +
        v10 * fx * (1 - fy) +
        v01 * (1 - fx) * fy +
        v11 * fx * fy
    )


# ============================================================================
# Spatial Hashing Helpers
# ============================================================================

@njit(cache=True)
def distance_squared(x1: float, y1: float, x2: float, y2: float) -> float:
    """Compute squared Euclidean distance."""
    dx = x1 - x2
    dy = y1 - y2
    return dx * dx + dy * dy


# ============================================================================
# Dispatcher Functions
# ============================================================================

def get_thomas_solver(n_rows: int):
    """Get the appropriate Thomas solver based on problem size."""
    if not HAS_NUMBA:
        # Fall back to pure numpy implementation
        from .environment import _thomas_solve_batch
        return _thomas_solve_batch

    # Use parallel version for larger grids
    if n_rows >= 50:
        return thomas_solve_batch_parallel
    else:
        return thomas_solve_batch_numba


# Pre-compile kernels on import (if numba available)
if HAS_NUMBA:
    # Trigger compilation with dummy data
    _dummy_a = np.zeros((2, 2), dtype=np.float64)
    _dummy_b = np.zeros((2, 3), dtype=np.float64)
    _dummy_c = np.zeros((2, 2), dtype=np.float64)
    _dummy_d = np.zeros((2, 3), dtype=np.float64)
    try:
        _ = thomas_solve_batch_numba(_dummy_a, _dummy_b, _dummy_c, _dummy_d)
    except Exception:
        pass  # Compilation may fail with invalid shapes, that's OK
