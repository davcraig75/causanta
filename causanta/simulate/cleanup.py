"""Cleanup utilities for CAUSANTA simulation outputs.

Reduces disk usage by removing intermediate files while preserving
key snapshots for analysis and visualization.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Callable


def cleanup_simulation_outputs(
    output_dir: Path,
    keep_fraction: float = 0.01,
    keep_first: bool = True,
    keep_last: bool = True,
    dry_run: bool = False,
    verbose: bool = True,
) -> dict[str, int]:
    """Remove intermediate TSV files, keeping only a fraction for analysis.

    This significantly reduces disk usage while preserving enough data for:
    - Time-series analysis
    - Animation generation
    - Final state analysis

    Args:
        output_dir: Directory containing simulation outputs
        keep_fraction: Fraction of files to keep (0.01 = 1%, 0.1 = 10%)
        keep_first: Always keep the first timestep (t=0)
        keep_last: Always keep the last timestep
        dry_run: If True, only report what would be deleted without deleting
        verbose: Print progress messages

    Returns:
        Dict with counts: {'cells_kept', 'cells_deleted', 'env_kept', 'env_deleted'}

    Example:
        >>> cleanup_simulation_outputs(Path("output/run_20260413"), keep_fraction=0.02)
        {'cells_kept': 20, 'cells_deleted': 980, 'env_kept': 20, 'env_deleted': 980}
    """
    stats = {
        'cells_kept': 0,
        'cells_deleted': 0,
        'env_kept': 0,
        'env_deleted': 0,
        'bytes_freed': 0,
    }

    # Find all timestamped files
    cells_files = sorted(output_dir.glob('cells_t*.tsv'))
    env_files = sorted(output_dir.glob('environment_t*.tsv'))

    if not cells_files:
        if verbose:
            print(f"No cells_t*.tsv files found in {output_dir}")
        return stats

    # Extract timesteps and determine which to keep
    def extract_time(path: Path) -> int:
        match = re.search(r't(\d+)', path.stem)
        return int(match.group(1)) if match else 0

    timesteps = sorted(set(extract_time(f) for f in cells_files))
    n_total = len(timesteps)
    n_keep = max(1, int(n_total * keep_fraction))

    # Calculate step size for uniform sampling
    if n_keep >= n_total:
        # Keep all
        keep_indices = set(range(n_total))
    else:
        step = n_total / n_keep
        keep_indices = set(int(i * step) for i in range(n_keep))

    # Always keep first and last if requested
    if keep_first:
        keep_indices.add(0)
    if keep_last:
        keep_indices.add(n_total - 1)

    keep_timesteps = set(timesteps[i] for i in keep_indices)

    if verbose:
        print(f"Cleanup: keeping {len(keep_timesteps)}/{n_total} timesteps ({len(keep_timesteps)/n_total*100:.1f}%)")
        print(f"  Keeping: t=0, ..., t={max(keep_timesteps)}")

    # Process cells files
    for f in cells_files:
        t = extract_time(f)
        if t in keep_timesteps:
            stats['cells_kept'] += 1
        else:
            stats['bytes_freed'] += f.stat().st_size
            if not dry_run:
                f.unlink()
            stats['cells_deleted'] += 1

    # Process environment files
    for f in env_files:
        t = extract_time(f)
        if t in keep_timesteps:
            stats['env_kept'] += 1
        else:
            stats['bytes_freed'] += f.stat().st_size
            if not dry_run:
                f.unlink()
            stats['env_deleted'] += 1

    if verbose:
        action = "Would delete" if dry_run else "Deleted"
        print(f"  {action}: {stats['cells_deleted']} cells files, {stats['env_deleted']} environment files")
        print(f"  Space freed: {stats['bytes_freed'] / 1024 / 1024:.1f} MB")

    return stats


def organize_output_directory(
    source_dir: Path,
    create_subdirs: bool = True,
) -> Path:
    """Organize simulation outputs into a clean directory structure.

    Creates subdirectories:
    - data/     - TSV files, logs, params
    - figures/  - Generated plots
    - reports/  - HTML reports
    - animations/ - Movie files

    Args:
        source_dir: Directory containing simulation outputs
        create_subdirs: Create subdirectories if they don't exist

    Returns:
        Path to the organized directory
    """
    if create_subdirs:
        (source_dir / "data").mkdir(exist_ok=True)
        (source_dir / "figures").mkdir(exist_ok=True)
        (source_dir / "reports").mkdir(exist_ok=True)
        (source_dir / "animations").mkdir(exist_ok=True)

    # Move files to appropriate subdirectories
    for f in source_dir.glob("*.tsv"):
        dest = source_dir / "data" / f.name
        if not dest.exists():
            f.rename(dest)

    for f in source_dir.glob("*.log"):
        dest = source_dir / "data" / f.name
        if not dest.exists():
            f.rename(dest)

    for f in source_dir.glob("params.json"):
        dest = source_dir / "data" / f.name
        if not dest.exists():
            f.rename(dest)

    for f in source_dir.glob("*.png"):
        dest = source_dir / "figures" / f.name
        if not dest.exists():
            f.rename(dest)

    for f in source_dir.glob("*.html"):
        if "movie" in f.name.lower() or "animation" in f.name.lower() or "growth" in f.name.lower():
            dest = source_dir / "animations" / f.name
        else:
            dest = source_dir / "reports" / f.name
        if not dest.exists():
            f.rename(dest)

    for f in source_dir.glob("*.vl.json"):
        dest = source_dir / "figures" / f.name
        if not dest.exists():
            f.rename(dest)

    return source_dir


def get_output_summary(output_dir: Path) -> dict:
    """Get summary statistics of output directory contents.

    Returns:
        Dict with file counts and sizes by category
    """
    summary = {
        'total_files': 0,
        'total_size_mb': 0,
        'categories': {},
    }

    categories = {
        'cells_tsv': list(output_dir.rglob('cells_t*.tsv')),
        'env_tsv': list(output_dir.rglob('environment_t*.tsv')),
        'lineage': list(output_dir.rglob('lineage*.tsv')),
        'html': list(output_dir.rglob('*.html')),
        'json': list(output_dir.rglob('*.json')),
        'log': list(output_dir.rglob('*.log')),
        'png': list(output_dir.rglob('*.png')),
    }

    for cat_name, files in categories.items():
        size = sum(f.stat().st_size for f in files if f.exists())
        summary['categories'][cat_name] = {
            'count': len(files),
            'size_mb': size / 1024 / 1024,
        }
        summary['total_files'] += len(files)
        summary['total_size_mb'] += size / 1024 / 1024

    return summary


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="Cleanup CAUSANTA simulation outputs")
    parser.add_argument("output_dir", type=Path, help="Output directory to clean")
    parser.add_argument("--keep", type=float, default=0.01, help="Fraction to keep (default: 0.01)")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be deleted")
    parser.add_argument("--organize", action="store_true", help="Organize into subdirectories")

    args = parser.parse_args()

    if args.organize:
        organize_output_directory(args.output_dir)
        print(f"Organized: {args.output_dir}")

    stats = cleanup_simulation_outputs(
        args.output_dir,
        keep_fraction=args.keep,
        dry_run=args.dry_run,
    )

    summary = get_output_summary(args.output_dir)
    print(f"\nRemaining: {summary['total_files']} files, {summary['total_size_mb']:.1f} MB")
