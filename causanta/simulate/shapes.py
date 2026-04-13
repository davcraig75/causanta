"""SVG shape geometry utilities for collision detection.

Parses SVG path strings and computes bounding boxes, effective radii,
and performs shape-aware collision detection.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from functools import lru_cache
from typing import Iterator

import numpy as np


@dataclass(frozen=True)
class ShapeBounds:
    """Axis-aligned bounding box for a shape in normalized [-1, 1] coordinates."""

    min_x: float
    max_x: float
    min_y: float
    max_y: float

    @property
    def width(self) -> float:
        return self.max_x - self.min_x

    @property
    def height(self) -> float:
        return self.max_y - self.min_y

    @property
    def center_x(self) -> float:
        return (self.min_x + self.max_x) / 2

    @property
    def center_y(self) -> float:
        return (self.min_y + self.max_y) / 2

    @property
    def max_radius(self) -> float:
        """Maximum radius from center to any corner."""
        corners = [
            (self.min_x - self.center_x, self.min_y - self.center_y),
            (self.max_x - self.center_x, self.min_y - self.center_y),
            (self.min_x - self.center_x, self.max_y - self.center_y),
            (self.max_x - self.center_x, self.max_y - self.center_y),
        ]
        return max(math.sqrt(cx * cx + cy * cy) for cx, cy in corners)


def _tokenize_path(path: str) -> Iterator[str]:
    """Tokenize SVG path into commands and numbers."""
    pattern = r"([MmLlHhVvCcSsQqTtAaZz])|(-?\d*\.?\d+(?:[eE][+-]?\d+)?)"
    for match in re.finditer(pattern, path):
        if match.group():
            yield match.group()


def _parse_numbers(tokens: list[str], start: int, count: int) -> tuple[list[float], int]:
    """Parse count numbers from tokens starting at index."""
    nums = []
    i = start
    while len(nums) < count and i < len(tokens):
        try:
            nums.append(float(tokens[i]))
            i += 1
        except ValueError:
            break
    return nums, i


@lru_cache(maxsize=64)
def parse_svg_path_bounds(path: str) -> ShapeBounds:
    """Parse SVG path and compute bounding box.

    Handles M, L, H, V, C, S, Q, T, A, Z commands (absolute only for now).
    Returns bounds in normalized coordinates.
    """
    if not path or not path.strip():
        return ShapeBounds(-0.5, 0.5, -0.5, 0.5)

    tokens = list(_tokenize_path(path))
    if not tokens:
        return ShapeBounds(-0.5, 0.5, -0.5, 0.5)

    points: list[tuple[float, float]] = []
    x, y = 0.0, 0.0
    start_x, start_y = 0.0, 0.0
    i = 0

    while i < len(tokens):
        cmd = tokens[i]
        if not cmd.isalpha():
            i += 1
            continue

        i += 1

        if cmd == "M":  # Move to
            nums, i = _parse_numbers(tokens, i, 2)
            if len(nums) >= 2:
                x, y = nums[0], nums[1]
                start_x, start_y = x, y
                points.append((x, y))

        elif cmd == "m":  # Relative move
            nums, i = _parse_numbers(tokens, i, 2)
            if len(nums) >= 2:
                x += nums[0]
                y += nums[1]
                start_x, start_y = x, y
                points.append((x, y))

        elif cmd == "L":  # Line to
            nums, i = _parse_numbers(tokens, i, 2)
            while len(nums) >= 2:
                x, y = nums[0], nums[1]
                points.append((x, y))
                nums, i = _parse_numbers(tokens, i, 2)

        elif cmd == "l":  # Relative line
            nums, i = _parse_numbers(tokens, i, 2)
            while len(nums) >= 2:
                x += nums[0]
                y += nums[1]
                points.append((x, y))
                nums, i = _parse_numbers(tokens, i, 2)

        elif cmd == "H":  # Horizontal line
            nums, i = _parse_numbers(tokens, i, 1)
            if len(nums) >= 1:
                x = nums[0]
                points.append((x, y))

        elif cmd == "h":  # Relative horizontal
            nums, i = _parse_numbers(tokens, i, 1)
            if len(nums) >= 1:
                x += nums[0]
                points.append((x, y))

        elif cmd == "V":  # Vertical line
            nums, i = _parse_numbers(tokens, i, 1)
            if len(nums) >= 1:
                y = nums[0]
                points.append((x, y))

        elif cmd == "v":  # Relative vertical
            nums, i = _parse_numbers(tokens, i, 1)
            if len(nums) >= 1:
                y += nums[0]
                points.append((x, y))

        elif cmd == "C":  # Cubic bezier
            nums, i = _parse_numbers(tokens, i, 6)
            while len(nums) >= 6:
                # Sample control points and endpoint
                points.append((nums[0], nums[1]))
                points.append((nums[2], nums[3]))
                x, y = nums[4], nums[5]
                points.append((x, y))
                nums, i = _parse_numbers(tokens, i, 6)

        elif cmd == "c":  # Relative cubic
            nums, i = _parse_numbers(tokens, i, 6)
            while len(nums) >= 6:
                points.append((x + nums[0], y + nums[1]))
                points.append((x + nums[2], y + nums[3]))
                x += nums[4]
                y += nums[5]
                points.append((x, y))
                nums, i = _parse_numbers(tokens, i, 6)

        elif cmd == "Q":  # Quadratic bezier
            nums, i = _parse_numbers(tokens, i, 4)
            while len(nums) >= 4:
                points.append((nums[0], nums[1]))
                x, y = nums[2], nums[3]
                points.append((x, y))
                nums, i = _parse_numbers(tokens, i, 4)

        elif cmd == "q":  # Relative quadratic
            nums, i = _parse_numbers(tokens, i, 4)
            while len(nums) >= 4:
                points.append((x + nums[0], y + nums[1]))
                x += nums[2]
                y += nums[3]
                points.append((x, y))
                nums, i = _parse_numbers(tokens, i, 4)

        elif cmd == "A":  # Arc
            nums, i = _parse_numbers(tokens, i, 7)
            while len(nums) >= 7:
                # rx, ry are radii; add points along arc envelope
                rx, ry = abs(nums[0]), abs(nums[1])
                ex, ey = nums[5], nums[6]
                # Include arc envelope (conservative)
                points.append((x - rx, y - ry))
                points.append((x + rx, y + ry))
                points.append((ex - rx, ey - ry))
                points.append((ex + rx, ey + ry))
                x, y = ex, ey
                points.append((x, y))
                nums, i = _parse_numbers(tokens, i, 7)

        elif cmd == "a":  # Relative arc
            nums, i = _parse_numbers(tokens, i, 7)
            while len(nums) >= 7:
                rx, ry = abs(nums[0]), abs(nums[1])
                ex, ey = x + nums[5], y + nums[6]
                points.append((x - rx, y - ry))
                points.append((x + rx, y + ry))
                points.append((ex - rx, ey - ry))
                points.append((ex + rx, ey + ry))
                x, y = ex, ey
                points.append((x, y))
                nums, i = _parse_numbers(tokens, i, 7)

        elif cmd in ("Z", "z"):  # Close path
            x, y = start_x, start_y
            points.append((x, y))

        elif cmd in ("S", "s", "T", "t"):
            # Smooth curves - parse appropriate number of coordinates
            count = 4 if cmd.upper() == "S" else 2
            nums, i = _parse_numbers(tokens, i, count)
            if cmd.islower():
                while len(nums) >= count:
                    for j in range(0, count, 2):
                        if j + 1 < count:
                            if cmd.islower():
                                points.append((x + nums[j], y + nums[j + 1]))
                            else:
                                points.append((nums[j], nums[j + 1]))
                    if count >= 2:
                        if cmd.islower():
                            x += nums[count - 2]
                            y += nums[count - 1]
                        else:
                            x, y = nums[count - 2], nums[count - 1]
                    nums, i = _parse_numbers(tokens, i, count)
            else:
                while len(nums) >= count:
                    for j in range(0, count, 2):
                        if j + 1 < count:
                            points.append((nums[j], nums[j + 1]))
                    x, y = nums[count - 2], nums[count - 1]
                    nums, i = _parse_numbers(tokens, i, count)

    if not points:
        return ShapeBounds(-0.5, 0.5, -0.5, 0.5)

    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    return ShapeBounds(min(xs), max(xs), min(ys), max(ys))


def get_effective_radius(shape_path: str, diameter: float) -> float:
    """Get effective collision radius for a shape.

    The shape_path is in normalized coordinates [-1, 1].
    Returns radius in micrometers scaled by diameter.
    """
    bounds = parse_svg_path_bounds(shape_path)
    # Max extent from origin determines collision radius
    # Multiply by half the diameter to convert to actual size
    return bounds.max_radius * (diameter / 2.0)


def get_oriented_bounds(
    shape_path: str,
    diameter: float,
    angle: float,
) -> tuple[float, float, float, float]:
    """Get axis-aligned bounds for a rotated shape.

    Returns (half_width, half_height) of the rotated bounding box.
    """
    bounds = parse_svg_path_bounds(shape_path)
    scale = diameter / 2.0

    # Corner points of the bounding box
    corners = [
        (bounds.min_x * scale, bounds.min_y * scale),
        (bounds.max_x * scale, bounds.min_y * scale),
        (bounds.min_x * scale, bounds.max_y * scale),
        (bounds.max_x * scale, bounds.max_y * scale),
    ]

    # Rotate corners
    cos_a, sin_a = math.cos(angle), math.sin(angle)
    rotated_x = []
    rotated_y = []
    for cx, cy in corners:
        rx = cx * cos_a - cy * sin_a
        ry = cx * sin_a + cy * cos_a
        rotated_x.append(rx)
        rotated_y.append(ry)

    half_w = max(abs(x) for x in rotated_x)
    half_h = max(abs(y) for y in rotated_y)
    return half_w, half_h, min(rotated_x), min(rotated_y)


def check_shape_collision(
    x1: float,
    y1: float,
    shape1: str,
    diameter1: float,
    angle1: float,
    x2: float,
    y2: float,
    shape2: str,
    diameter2: float,
    angle2: float,
    margin: float = 0.0,
) -> bool:
    """Check if two shapes overlap using axis-aligned bounding box check.

    Uses a two-phase approach:
    1. Quick distance check with average radii (filters obvious non-collisions)
    2. AABB intersection test for closer shapes

    Returns True if shapes likely overlap.
    """
    # Get rotated half-extents for each shape
    hw1, hh1, _, _ = get_oriented_bounds(shape1, diameter1, angle1)
    hw2, hh2, _, _ = get_oriented_bounds(shape2, diameter2, angle2)

    # Add margin
    hw1 += margin / 2
    hh1 += margin / 2
    hw2 += margin / 2
    hh2 += margin / 2

    # Compute center distance
    dx = abs(x2 - x1)
    dy = abs(y2 - y1)

    # AABB intersection: shapes collide if centers are closer than sum of half-extents
    # This is the correct AABB test: check each axis independently
    x_overlap = dx < (hw1 + hw2)
    y_overlap = dy < (hh1 + hh2)

    return x_overlap and y_overlap


def check_circle_shape_collision(
    x1: float,
    y1: float,
    radius1: float,
    x2: float,
    y2: float,
    shape2: str,
    diameter2: float,
    angle2: float,
    margin: float = 0.0,
) -> bool:
    """Check collision between a circle and a shape."""
    hw2, hh2, _, _ = get_oriented_bounds(shape2, diameter2, angle2)

    # Add margin
    r1 = radius1 + margin / 2
    hw2 += margin / 2
    hh2 += margin / 2

    # For circle vs AABB: check if circle center is within expanded AABB
    dx = abs(x2 - x1)
    dy = abs(y2 - y1)

    return dx < (r1 + hw2) and dy < (r1 + hh2)


# Pre-compute bounds for default shapes
DEFAULT_SHAPE_BOUNDS: dict[str, ShapeBounds] = {}


def get_cached_bounds(shape_path: str) -> ShapeBounds:
    """Get bounds from cache or compute."""
    if shape_path in DEFAULT_SHAPE_BOUNDS:
        return DEFAULT_SHAPE_BOUNDS[shape_path]
    return parse_svg_path_bounds(shape_path)


def transform_svg_path(
    path: str,
    scale: float,
    translate_x: float,
    translate_y: float,
    rotate: float = 0.0,
) -> str:
    """Transform an SVG path by scale, translation, and rotation.

    Used for rendering cells at their actual positions.
    """
    if not path:
        return path

    # For efficiency, return a transform string rather than rewriting path
    # This can be used with SVG transform attribute
    parts = []
    if translate_x != 0 or translate_y != 0:
        parts.append(f"translate({translate_x:.1f},{translate_y:.1f})")
    if rotate != 0:
        parts.append(f"rotate({math.degrees(rotate):.1f})")
    if scale != 1:
        parts.append(f"scale({scale:.2f})")

    return " ".join(parts) if parts else ""
