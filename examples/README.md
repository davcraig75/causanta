# CAUSANTA Examples

This directory contains example scripts demonstrating how to use CAUSANTA for tumor simulation and causal inference.

## Scripts

### `analyze_simulation.py`

Demonstrates how to analyze simulation output using the causal inference tools:

```bash
python examples/analyze_simulation.py output/run_20260413_093442/data/
```

Features:
- Load simulation output (cells, lineage, environment)
- Run instrumental variable regression using ecDNA as instrument
- Estimate causal effects of EGFR expression on migration
- Compare estimated vs. true (configured) causal parameters

## Quick Start

1. **Run a simulation:**
   ```bash
   python -m causanta.simulate.core --hours 200
   ```

2. **Analyze the output:**
   ```bash
   python examples/analyze_simulation.py output/run_*/data/
   ```

## See Also

- [Tutorial: Causal Inference with CAUSANTA](../docs/tutorial.html)
- [Analysis Guide](../docs/analysis_guide.html)
- [Immune System Documentation](../docs/immune_system.html)
