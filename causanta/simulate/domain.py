"""Spatial domain, coordinate system, and grid operations for CAUSANTA.

Performance: Uses Numba JIT compilation when available for hot paths.
"""

from __future__ import annotations

import numpy as np

from .config import DomainConfig

# Import Numba-accelerated kernels if available
try:
    from ._numba_kernels import (
        HAS_NUMBA,
        bilinear_interpolate as _numba_bilinear,
        compute_gradient_at_point as _numba_gradient,
    )
except ImportError:
    HAS_NUMBA = False
    _numba_bilinear = None
    _numba_gradient = None


class Domain:
    """2D tissue domain with environment grid and coordinate transforms."""

    def __init__(self, config: DomainConfig):
        self.width_um = config.width_um
        self.height_um = config.height_um
        self.env_grid_um = config.env_grid_um
        self.nx = config.width_um // config.env_grid_um
        self.ny = config.height_um // config.env_grid_um

    def um_to_grid(self, x_um: int, y_um: int) -> tuple[int, int]:
        """Convert micrometer position to grid indices (col, row).

        Returns clamped indices within [0, nx-1] x [0, ny-1].
        """
        gi = min(max(int(x_um) // self.env_grid_um, 0), self.nx - 1)
        gj = min(max(int(y_um) // self.env_grid_um, 0), self.ny - 1)
        return gi, gj

    def grid_to_um(self, gi: int, gj: int) -> tuple[float, float]:
        """Convert grid indices to micrometer center position."""
        x = (gi + 0.5) * self.env_grid_um
        y = (gj + 0.5) * self.env_grid_um
        return x, y

    def bilinear_interpolate(self, field: np.ndarray, x_um: float, y_um: float) -> float:
        """Sample a 2D field at a continuous position using bilinear interpolation.

        field has shape (ny, nx) with field[row, col] indexing.
        """
        # Use Numba version if available
        if HAS_NUMBA and _numba_bilinear is not None:
            return _numba_bilinear(field, x_um, y_um, float(self.env_grid_um))

        # Pure Python fallback
        return self._bilinear_interpolate_python(field, x_um, y_um)

    def _bilinear_interpolate_python(self, field: np.ndarray, x_um: float, y_um: float) -> float:
        """Pure Python bilinear interpolation (fallback)."""
        # Convert to continuous grid coordinates (center of voxel 0 is at 0.5*grid_um)
        gx = x_um / self.env_grid_um - 0.5
        gy = y_um / self.env_grid_um - 0.5

        # Floor indices
        x0 = int(np.floor(gx))
        y0 = int(np.floor(gy))

        # Fractional parts
        fx = gx - x0
        fy = gy - y0

        # Clamp to valid range
        x0c = min(max(x0, 0), self.nx - 2)
        y0c = min(max(y0, 0), self.ny - 2)
        x1c = x0c + 1
        y1c = y0c + 1

        # Bilinear interpolation: field[row, col]
        v00 = field[y0c, x0c]
        v10 = field[y0c, x1c]
        v01 = field[y1c, x0c]
        v11 = field[y1c, x1c]

        return float(
            v00 * (1 - fx) * (1 - fy)
            + v10 * fx * (1 - fy)
            + v01 * (1 - fx) * fy
            + v11 * fx * fy
        )

    def compute_gradient(self, field: np.ndarray, x_um: float, y_um: float) -> tuple[float, float]:
        """Compute gradient (dF/dx, dF/dy) at a position via central differences.

        Returns gradient in units of field_value / um.
        """
        # Use Numba version if available
        if HAS_NUMBA and _numba_gradient is not None:
            return _numba_gradient(field, x_um, y_um, float(self.env_grid_um))

        # Pure Python fallback
        return self._compute_gradient_python(field, x_um, y_um)

    def _compute_gradient_python(self, field: np.ndarray, x_um: float, y_um: float) -> tuple[float, float]:
        """Pure Python gradient computation (fallback)."""
        gi, gj = self.um_to_grid(int(x_um), int(y_um))

        # Central differences with boundary clamping
        gi_lo = max(gi - 1, 0)
        gi_hi = min(gi + 1, self.nx - 1)
        gj_lo = max(gj - 1, 0)
        gj_hi = min(gj + 1, self.ny - 1)

        dx_grid = gi_hi - gi_lo
        dy_grid = gj_hi - gj_lo

        dfdx = 0.0
        if dx_grid > 0:
            dfdx = (field[gj, gi_hi] - field[gj, gi_lo]) / (dx_grid * self.env_grid_um)

        dfdy = 0.0
        if dy_grid > 0:
            dfdy = (field[gj_hi, gi] - field[gj_lo, gi]) / (dy_grid * self.env_grid_um)

        return dfdx, dfdy

    def is_in_bounds(self, x_um: int, y_um: int) -> bool:
        """Check if position is within domain boundaries."""
        return 0 <= x_um < self.width_um and 0 <= y_um < self.height_um

    def clamp_position(self, x_um: int, y_um: int) -> tuple[int, int]:
        """Clamp position to domain boundaries."""
        return (
            min(max(x_um, 0), self.width_um - 1),
            min(max(y_um, 0), self.height_um - 1),
        )
