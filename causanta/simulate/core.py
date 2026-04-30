"""Main simulation loop orchestration for CAUSANTA.

Implements the 6-phase simulation loop from Section 8 of the framework:
1. Environment diffusion
2. Cell behaviors (random-shuffled order)
3. Angiogenesis
4. Immune recruitment
5. Cleanup
6. Output
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

from .angiogenesis import AngiogenesisManager
from .behaviors import (
    LineageRecord,
    advance_cell_cycle,
    apply_necrosis,
    check_apoptosis,
    check_immune_kill,
    check_necrosis,
    check_state_transitions,
    clear_necrotic_cells,
    compute_migration,
    sample_environment,
    update_effective_rates,
    update_immune_activation,
    update_immune_exhaustion,
)
from .cells import (
    ENDOTHELIAL,
    MICROGLIA,
    NECROTIC,
    RECRUITED_IMMUNE,
    TUMOR,
    CellPopulation,
    create_cell,
)
from .config import SimulationConfig, load_config
from .domain import Domain
from .environment import (
    DiffusionSolver,
    EnvironmentFields,
    accumulate_sources_and_sinks,
)
from .initialization import (
    equilibrate_environment,
    generate_vascular_network,
    populate_normal_tissue,
    seed_tumor,
)
from .cleanup import cleanup_data_keep_final_only, compress_movie_html
from .io import (
    copy_params_json,
    create_output_directory,
    write_cells_tsv,
    write_environment_tsv,
    write_lineage_tsv,
    write_summary_log,
)
from .movie import generate_movie_html
from .reporting import generate_report
from .viewer import generate_standalone_html, generate_vega_json
from .visualization import write_html_viewer, write_vega_spec


class Simulation:
    """CAUSANTA simulation engine."""

    def __init__(self, config_path: str | Path):
        self.config_path = Path(config_path)
        self.config = load_config(config_path)
        self.rng = np.random.default_rng(self.config.rng_seed)

        # Core components
        self.domain = Domain(self.config.domain)
        self.env = EnvironmentFields(self.domain, self.config.environment)
        self.solver = DiffusionSolver(self.domain, self.config.environment)
        self.population = CellPopulation(self.domain)
        self.angiogenesis = AngiogenesisManager(self.config)

        # State tracking
        self.current_time_hr: float = 0.0
        self.lineage_records: list[LineageRecord] = []
        self.output_dir: Path | None = None

    def initialize(self) -> None:
        """Initialize tissue, environment, and output directory (spec Section 8)."""
        print("CAUSANTA Simulation Engine v0.1.0")
        print(f"  Domain: {self.domain.width_um} x {self.domain.height_um} um")
        print(f"  Env grid: {self.domain.env_grid_um} um ({self.domain.nx} x {self.domain.ny})")
        print(f"  Total time: {self.config.time.total_hours} hours")
        print()

        # Create output directory
        if self.config.output_dir:
            # Use explicit output directory from config
            self.output_dir = Path(self.config.output_dir)
            self.output_dir.mkdir(parents=True, exist_ok=True)
            # Create subdirectories
            (self.output_dir / "data").mkdir(exist_ok=True)
            (self.output_dir / "figures").mkdir(exist_ok=True)
            (self.output_dir / "reports").mkdir(exist_ok=True)
            (self.output_dir / "animations").mkdir(exist_ok=True)
        else:
            self.output_dir = create_output_directory()
        copy_params_json(self.config_path, self.output_dir)
        print(f"  Output: {self.output_dir}")
        print()

        # Phase 1: Generate vascular network
        print("Initializing tissue...")
        generate_vascular_network(
            self.domain, self.population, self.env, self.config, self.rng
        )
        endo_count = sum(1 for _ in self.population.iter_by_type(ENDOTHELIAL))
        print(f"  Vascular network: {endo_count} endothelial cells")

        # Phase 2: Populate normal cells
        populate_normal_tissue(self.domain, self.population, self.config, self.rng)
        counts = self.population.count_by_type()
        total = self.population.count()
        print(f"  Total cells: {total}")
        for tid, cnt in sorted(counts.items()):
            name = self.config.cell_types.get(tid)
            tname = name.type_name if name else f"Type {tid}"
            print(f"    {tname}: {cnt}")

        # Phase 3: Equilibrate environment (BEFORE tumor seeding)
        # This ensures glucose/O2 reach steady state before tumor starts consuming
        equilibrate_environment(
            self.env, self.solver, self.population, self.domain, self.config
        )
        print(f"  Equilibrated: mean O2 = {self.env.O2.mean():.1f} mmHg")
        print()

        # Phase 4: Seed tumor (AFTER equilibration)
        seed_tumor(self.population, self.config, self.rng, 0.0)
        tumor_count = sum(1 for _ in self.population.iter_by_type(TUMOR))
        print(f"  Tumor cells seeded: {tumor_count}")
        print()

        # Write initial state
        self._write_output(0)
        self._write_lineage()

    def run(self) -> None:
        """Run the full simulation."""
        self.initialize()

        total = self.config.time.total_hours
        t_start = time.time()

        print(f"Running simulation for {total} hours...")
        for step in range(1, total + 1):
            self.current_time_hr = float(step)

            # Phase 1: Environment diffusion
            self._step_environment()

            # Phase 2: Cell behaviors
            self._step_cells()

            # Phase 3: Angiogenesis
            self._step_angiogenesis()

            # Phase 4: Immune recruitment
            self._step_immune_recruitment()

            # Phase 5: Cleanup
            self._step_cleanup()

            # Phase 6: Output
            if step % self.config.time.output_interval_hr == 0:
                self._write_output(step)

            # Write any new lineage records
            if self.lineage_records:
                self._write_lineage()

            # Progress report
            if step % max(1, total // 20) == 0 or step == total:
                elapsed = time.time() - t_start
                counts = self.population.count_by_type()
                total_cells = self.population.count()
                tumor_n = counts.get(TUMOR, 0)
                necrotic_n = counts.get(NECROTIC, 0)
                mean_o2 = float(self.env.O2.mean())
                rate = step / elapsed if elapsed > 0 else 0
                print(
                    f"  t={step:>4d}h | cells={total_cells:>5d} | "
                    f"tumor={tumor_n:>4d} | necrotic={necrotic_n:>3d} | "
                    f"O2={mean_o2:.1f} | {rate:.1f} hr/s"
                )

        elapsed = time.time() - t_start
        print(f"\nSimulation complete in {elapsed:.1f}s")
        print(f"  Output: {self.output_dir}")

        # Write visualization and report
        if self.output_dir:
            # Determine output subdirectories
            reports_dir = self.output_dir / "reports"
            figures_dir = self.output_dir / "figures"
            data_dir = self.output_dir / "data"

            # Use subdirectories if they exist, otherwise use root
            if not reports_dir.exists():
                reports_dir = self.output_dir
            if not figures_dir.exists():
                figures_dir = self.output_dir
            if not data_dir.exists():
                data_dir = self.output_dir

            write_vega_spec(
                data_dir,
                self.config,
                total,
                self.config.time.output_interval_hr,
            )

            # Generate Vega JSON spec with embedded data (must be first for index.html)
            vega_json_path = generate_vega_json(
                population=self.population,
                env=self.env,
                domain=self.domain,
                config=self.config,
                time_hr=float(total),
                output_path=figures_dir / "simulation.vl.json",
            )
            print(f"  Vega spec: {vega_json_path}")

            # Load the spec and embed it in index.html for standalone viewing
            with open(vega_json_path) as f:
                vega_spec = json.load(f)
            write_html_viewer(reports_dir, vega_spec=vega_spec)
            print(f"  Visualization: {reports_dir / 'index.html'}")

            # Generate standalone HTML viewer with embedded Vega-Lite
            standalone_path = generate_standalone_html(
                population=self.population,
                env=self.env,
                domain=self.domain,
                config=self.config,
                time_hr=float(total),
                output_path=reports_dir / "viewer.html",
                title=f"CAUSANTA Simulation (t={total}h)",
            )
            print(f"  Standalone viewer: {standalone_path}")

            # Generate animated movie from time-series data
            print("Generating animation...")
            try:
                movie_path = generate_movie_html(
                    output_dir=data_dir,
                    movie_path=reports_dir / "movie.html",
                    domain_width=self.domain.width_um,
                    domain_height=self.domain.height_um,
                    fps=4,
                    max_frames=200,
                    render_mode="tumor_focus",
                )
                print(f"  Animation: {movie_path}")
            except ValueError as e:
                print(f"  Animation skipped: {e}")

            # Generate comprehensive report
            print("Generating report...")
            report_path = generate_report(reports_dir, self.config, elapsed, data_dir=data_dir)
            print(f"  Report: {report_path}")

    def _step_environment(self) -> None:
        """Phase 1: Diffusion sub-stepping for all substrates."""
        dt = self.config.time.dt_diffusion_hr
        n_substeps = max(1, int(1.0 / dt))

        sources, sinks = accumulate_sources_and_sinks(
            self.population, self.env, self.domain, self.config
        )
        self.solver.step(self.env, sources, sinks, dt, n_substeps)

    def _step_cells(self) -> None:
        """Phase 2: Process cells in random-shuffled order."""
        cells = self.population.shuffled_living(self.rng)
        hypoxia_thresh = self.config.environment.hypoxia_threshold_mmHg
        kill_radius = self.config.immune_recruitment.kill_radius_um

        to_remove: list[int] = []

        for cell in cells:
            if not cell.is_alive:
                continue

            tp = self.config.cell_types.get(cell.cell_type)
            if tp is None:
                continue

            # a. Sample local environment
            sample_environment(cell, self.env, self.domain, hypoxia_thresh)
            update_effective_rates(
                cell, tp, self.rng,
                egfr_hypoxia_upregulation=self.config.environment.egfr_hypoxia_upregulation,
            )

            # b. Death checks
            # Necrosis
            if check_necrosis(cell, tp):
                apply_necrosis(cell, self.config)
                continue

            # Apoptosis
            if check_apoptosis(cell, tp, self.rng):
                cell.is_alive = False
                to_remove.append(cell.cell_id)
                continue

            # c. State transitions
            check_state_transitions(
                cell, tp, self.env, self.domain, self.population, self.config, self.rng
            )

            # d. Proliferation
            if tp.can_divide:
                record = advance_cell_cycle(
                    cell, tp, self.population, self.config, self.rng, self.current_time_hr
                )
                if record:
                    self.lineage_records.append(record)

            # e. Migration
            compute_migration(
                cell, tp, self.env, self.domain, self.population, self.rng
            )

        # Immune cell state updates and killing
        # Process recruited immune cells
        for cell in self.population.iter_by_type(RECRUITED_IMMUNE):
            if not cell.is_alive:
                continue

            # Update activation based on tumor proximity
            update_immune_activation(cell, self.population, self.config)

            # Update exhaustion recovery
            update_immune_exhaustion(cell, self.config)

            # Attempt to kill nearby tumor cells
            nearby = self.population.get_neighbors(cell.x, cell.y, kill_radius)
            for target in nearby:
                if target.cell_type == TUMOR and target.is_alive:
                    if check_immune_kill(target, cell, self.config, self.rng):
                        target.is_alive = False
                        to_remove.append(target.cell_id)
                        break  # One kill attempt per step

        # Also check activated microglia (resident immune cells)
        for cell in self.population.iter_by_type(MICROGLIA):
            if not cell.is_alive:
                continue

            # Update activation state
            update_immune_activation(cell, self.population, self.config)

            # Only attempt kill if reactive/activated
            if not cell.is_reactive:
                continue

            # Update exhaustion recovery
            update_immune_exhaustion(cell, self.config)

            nearby = self.population.get_neighbors(cell.x, cell.y, kill_radius)
            for target in nearby:
                if target.cell_type == TUMOR and target.is_alive:
                    if check_immune_kill(target, cell, self.config, self.rng):
                        target.is_alive = False
                        to_remove.append(target.cell_id)
                        break  # One kill attempt per step

        # Remove dead cells
        for cid in to_remove:
            self.population.remove_cell(cid)

    def _step_angiogenesis(self) -> None:
        """Phase 3: Vascular sprouting."""
        self.angiogenesis.process(
            self.population, self.env, self.domain, self.rng, self.current_time_hr
        )

    def _step_immune_recruitment(self) -> None:
        """Phase 4: Spawn immune cells from vasculature."""
        recruit_cfg = self.config.immune_recruitment
        immune_params = self.config.cell_types.get(RECRUITED_IMMUNE)
        if immune_params is None:
            return

        # For each vascular position, check if chemokine-like signal is present
        # Approximate chemokine as VEGF + lactate proxy near tumor
        for gj in range(self.domain.ny):
            for gi in range(self.domain.nx):
                if self.env.vascular_density[gj, gi] < 0.3:
                    continue

                # Chemokine proxy: VEGF > 0 indicates tumor nearby
                local_vegf = self.env.VEGF[gj, gi]
                if local_vegf < recruit_cfg.chemokine_threshold:
                    continue

                if self.rng.random() > recruit_cfg.recruitment_rate_per_hr:
                    continue

                # Spawn immune cell at vascular position
                x_um, y_um = self.domain.grid_to_um(gi, gj)
                x, y = int(x_um), int(y_um)
                if not self.population.has_neighbor_within(
                    x, y, immune_params.diameter_mean_um * 0.6
                ):
                    cid = self.population.allocate_id()
                    cell = create_cell(
                        cell_id=cid,
                        parent_id=-1,
                        cell_type=RECRUITED_IMMUNE,
                        x=x,
                        y=y,
                        time_born_hr=self.current_time_hr,
                        type_params=immune_params,
                        rng=self.rng,
                        cell_cycle_phase="G0",
                    )
                    cell.activation_level = 0.5
                    self.population.add_cell(cell)

    def _step_cleanup(self) -> None:
        """Phase 5: Lysis of necrotic cells and update occupancy."""
        clear_necrotic_cells(self.population, self.config, self.rng)
        self._update_occupancy()

    def _update_occupancy(self) -> None:
        """Update environment grid occupancy from current cell positions."""
        self.env.clear_occupancy()
        for cell in self.population.iter_living():
            gi, gj = self.domain.um_to_grid(cell.x, cell.y)
            self.env.is_occupied[gj, gi] = True
            self.env.occupant_cell_id[gj, gi] = cell.cell_id
            self.env.occupant_type[gj, gi] = cell.cell_type

    def _write_output(self, step: int) -> None:
        """Phase 6: Write TSV files and summary statistics."""
        if self.output_dir is None:
            return

        # Determine data directory (use subdirectory if it exists)
        data_dir = self.output_dir / "data"
        if not data_dir.exists():
            data_dir = self.output_dir

        # Cell state
        cells_path = data_dir / f"cells_t{step:06d}.tsv"
        write_cells_tsv(self.population, cells_path)

        # Environment state
        env_path = data_dir / f"environment_t{step:06d}.tsv"
        write_environment_tsv(self.env, self.domain, env_path)

        # Summary statistics
        counts = self.population.count_by_type()
        stats = {
            "total_cells": self.population.count(),
            "tumor_cells": counts.get(TUMOR, 0),
            "necrotic_cells": counts.get(NECROTIC, 0),
            "immune_cells": counts.get(RECRUITED_IMMUNE, 0),
            "mean_O2_mmHg": float(self.env.O2.mean()),
            "min_O2_mmHg": float(self.env.O2.min()),
            "mean_glucose_mM": float(self.env.glucose.mean()),
            "max_VEGF_nM": float(self.env.VEGF.max()),
            "mean_lactate_mM": float(self.env.lactate.mean()),
        }
        summary_path = data_dir / "summary.log"
        write_summary_log(stats, summary_path, step, append=(step > 0))

    def _write_lineage(self) -> None:
        """Write pending lineage records and clear the buffer."""
        if self.output_dir is None or not self.lineage_records:
            return
        # Use data subdirectory if it exists
        data_dir = self.output_dir / "data"
        dest_dir = data_dir if data_dir.exists() else self.output_dir
        lineage_path = dest_dir / "lineage.tsv"
        append = lineage_path.exists()
        write_lineage_tsv(self.lineage_records, lineage_path, append=append)
        self.lineage_records.clear()


def main() -> None:
    """CLI entry point for running CAUSANTA simulations."""
    parser = argparse.ArgumentParser(
        description="CAUSANTA: Causal Analysis Using Somatic And Neighborhood Tissue Architecture"
    )
    parser.add_argument(
        "config",
        nargs="?",
        default=None,
        help="Path to JSON configuration file (default: built-in params/default.json)",
    )
    parser.add_argument(
        "--hours",
        type=int,
        default=None,
        help="Override total simulation hours",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Override RNG seed",
    )
    parser.add_argument(
        "--keep-all-snapshots",
        action="store_true",
        help=(
            "Keep every per-hour cells_t*/environment_t* TSV. "
            "Default: only the final timestep is retained at end of run."
        ),
    )
    parser.add_argument(
        "--no-compress-movie",
        action="store_true",
        help=(
            "Skip end-of-run gzip of reports/movie.html. "
            "Default: movie.html is replaced by movie.html.gz."
        ),
    )
    args = parser.parse_args()

    # Determine config path
    if args.config:
        config_path = Path(args.config)
    else:
        config_path = Path(__file__).parent / "params" / "default.json"

    if not config_path.exists():
        print(f"Error: config file not found: {config_path}", file=sys.stderr)
        sys.exit(1)

    sim = Simulation(config_path)

    # Apply CLI overrides
    tmp_path = None
    if args.hours is not None:
        # Rebuild config with overridden hours — use mutable workaround
        import json
        import tempfile
        with open(config_path) as f:
            data = json.load(f)
        data.setdefault("time", {})["total_hours"] = args.hours
        if args.seed is not None:
            data["rng_seed"] = args.seed
        # Write temp config
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as tmp:
            json.dump(data, tmp)
            tmp_path = tmp.name
        sim = Simulation(tmp_path)
        sim.config_path = config_path  # Use original for copy
    elif args.seed is not None:
        sim.rng = np.random.default_rng(args.seed)

    try:
        sim.run()

        # Post-run space optimisation. Movie compression first (its input is the
        # untouched movie.html); snapshot cleanup second. Both run after every
        # downstream artifact (movie, viewer, report) has been generated, so
        # they only affect on-disk persistence, not the simulation results.
        if sim.output_dir is not None:
            reports_dir = sim.output_dir / "reports"
            if not reports_dir.exists():
                reports_dir = sim.output_dir
            data_dir = sim.output_dir / "data"
            if not data_dir.exists():
                data_dir = sim.output_dir

            if not args.no_compress_movie:
                compress_movie_html(reports_dir / "movie.html")

            if not args.keep_all_snapshots:
                cleanup_data_keep_final_only(data_dir)
    finally:
        # Clean up temp config file
        if tmp_path is not None:
            import os
            try:
                os.unlink(tmp_path)
            except OSError:
                pass


if __name__ == "__main__":
    main()
