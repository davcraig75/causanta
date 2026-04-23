"""Visualization module for CAUSANTA framework.

Provides Nature Methods-styled figures for publication.
"""

from .nature_style import (
    apply_nature_style,
    save_figure,
    NATURE_COLORS,
    NATURE_FIGSIZE,
)
from .causal_figures import (
    plot_dag,
    plot_segregation_validation,
    plot_first_stage,
    plot_iv_estimates,
    plot_power_curves,
    plot_heterogeneity,
    plot_sensitivity,
    create_main_figure_1,
    create_main_figure_2,
    create_main_figure_3,
    create_main_figure_4,
    create_main_figure_5,
    create_main_figure_6,
)

__all__ = [
    "apply_nature_style",
    "save_figure",
    "NATURE_COLORS",
    "NATURE_FIGSIZE",
    "plot_dag",
    "plot_segregation_validation",
    "plot_first_stage",
    "plot_iv_estimates",
    "plot_power_curves",
    "plot_heterogeneity",
    "plot_sensitivity",
    "create_main_figure_1",
    "create_main_figure_2",
    "create_main_figure_3",
    "create_main_figure_4",
    "create_main_figure_5",
    "create_main_figure_6",
]
