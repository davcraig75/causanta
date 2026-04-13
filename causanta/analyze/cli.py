"""Command-line interface for CAUSANTA causal analysis."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from .loader import load_simulation_output
from .report import generate_analysis_report


def main() -> None:
    """CLI entry point for running CAUSANTA causal analysis."""
    parser = argparse.ArgumentParser(
        description="CAUSANTA: Causal Analysis of Simulation Output"
    )
    parser.add_argument(
        "output_dir",
        help="Path to simulation output directory",
    )
    parser.add_argument(
        "--methods",
        type=str,
        default="all",
        help="Analysis methods to run: all, iv, discovery, matching, regression, sem, sensitivity",
    )
    parser.add_argument(
        "--dag-config",
        type=str,
        default=None,
        help="Path to DAG configuration JSON file",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Output path for analysis report (default: analysis_report.html in output_dir)",
    )
    parser.add_argument(
        "--format",
        type=str,
        choices=["html", "json", "both"],
        default="html",
        help="Output format for report",
    )

    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    if not output_dir.exists():
        print(f"Error: output directory not found: {output_dir}", file=sys.stderr)
        sys.exit(1)

    # Load data
    print(f"Loading simulation output from: {output_dir}")
    data = load_simulation_output(output_dir, load_all_timesteps=True)
    print(f"  Loaded {data.n_timesteps} timesteps")
    print(f"  {len(data.tumor_cells)} tumor cells at final timestep")
    print(f"  {data.n_divisions} division events")

    # Parse methods
    if args.methods == "all":
        methods = ["iv", "discovery", "matching", "regression", "sem", "sensitivity"]
    else:
        methods = [m.strip() for m in args.methods.split(",")]

    # Load DAG config if provided
    dag_config = None
    if args.dag_config:
        with open(args.dag_config) as f:
            dag_config = json.load(f)

    # Generate report
    output_path = args.output
    if output_path is None:
        output_path = output_dir / "analysis_report.html"
    else:
        output_path = Path(output_path)

    print(f"\nRunning analysis with methods: {', '.join(methods)}")
    report_path = generate_analysis_report(
        data=data,
        output_path=output_path,
        methods=methods,
        dag_config=dag_config,
        include_json=(args.format in ["json", "both"]),
    )

    print(f"\nAnalysis complete!")
    print(f"  Report: {report_path}")


if __name__ == "__main__":
    main()
