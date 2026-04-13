# CAUSANTA Codex Response

## Systematic Response to Deep Review (codex.md)

Date: 2026-04-13

Each issue is addressed with: diagnosis, action taken, and resolution status.

---

## Issue 1: IV benchmark invalid — analyzer uses ecDNA as both instrument AND treatment

**Diagnosis:** Confirmed. `iv.py:207` sets `D = Z` (both ecDNA_count). The loader also drops `egfr_expression` from the TSV parse. This makes 2SLS a regression of a variable on itself.

**Action:**
- `loader.py`: Added `egfr_expression` to the column parser
- `iv.py`: Changed treatment variable D to `egfr_expression` when available, falling back to ecDNA only with a warning
- Discovery default variables now include `egfr_expression`

**Status:** RESOLVED

---

## Issue 2: Simulator violates exclusion restriction — apoptosis and division use ecDNA directly

**Diagnosis:** Confirmed.
- `check_apoptosis()` passes `cell.ecDNA_count` to `modulate_apoptosis_rate()` instead of `cell.egfr_expression`
- G0→G1 transition passes `cell.ecDNA_count` to `modulate_division_time()` instead of `cell.egfr_expression`
- `_execute_division()` correctly uses `parent.egfr_expression` after calling `update_effective_rates()`

**Action:**
- `check_apoptosis()`: Changed to use `cell.egfr_expression`
- G0→G1 `advance_cell_cycle()`: Changed to use `cell.egfr_expression`
- Both now route through the same causal chain: ecDNA → EGFR → phenotype

**Status:** RESOLVED

---

## Issue 3: Segregation analysis math inconsistent with replication model

**Diagnosis:** Confirmed. The simulator replicates ecDNA during S-phase (nearly doubling count) then segregates the replicated pool. So `daughter_ecDNA / parent_ecDNA_before` has expected value ≈ 0.95 (replication fidelity), not 0.5.

**Action:**
- `effects.py`: Segregation analysis now computes daughter fraction relative to `parent_ecDNA_before + parent_ecDNA_after` (the total post-replication pool), where expected fraction is 0.5
- Updated docstring to explain the replication-then-segregation model
- Tutorial and analysis guide documentation updated to reflect this

**Status:** RESOLVED

---

## Issue 4: Environment physics are unit-inconsistent

**Diagnosis:** Confirmed. The reaction terms mix mmHg, mM, amol/hr, and um² without conversion factors. The fields are phenomenological, not physically dimensional.

**Action:** This is a design limitation, not a bug to fix in code. Added explicit documentation:
- Comment block in `environment.py` source stating fields are phenomenological
- Updated analysis_guide.md with a "Caveats" section noting unit inconsistency
- This is acceptable for a benchmark prototype; calibrated units are a Priority 4 item

**Status:** RESOLVED (documented as known limitation)

---

## Issue 5: GES algorithm can produce cyclic graphs

**Diagnosis:** Confirmed. The forward phase adds edges greedily without checking for cycles. Can produce A→B and B→A simultaneously.

**Action:**
- Added `_would_create_cycle()` helper using DFS reachability check
- GES forward phase now skips any edge addition that would create a cycle
- Backward phase is unaffected (removing edges cannot create cycles)

**Status:** RESOLVED

---

## Issue 6: Immune chemotaxis not wired; max_kills_before_exhaustion unused

**Diagnosis:** Confirmed.
- `compute_immune_chemotaxis()` exists but `compute_migration()` never calls it
- `max_kills_before_exhaustion` is in config but not enforced in `check_immune_kill()`

**Action:**
- `compute_migration()`: Added call to `compute_immune_chemotaxis()` for immune cell types (RecruitedImmune and reactive Microglia)
- `check_immune_kill()`: Added hard cap — returns False when `kills_performed >= max_kills_before_exhaustion`

**Status:** RESOLVED

---

## Issue 7: Reporting/viewer output paths inconsistent

**Diagnosis:** Confirmed.
- `index.html` references `visualization.vg.json` but code writes `simulation.vl.json`
- Report lands in `data/` instead of `reports/`

**Action:**
- `visualization.py`: Changed the vegaEmbed reference to `../figures/simulation.vl.json` (correct relative path from reports/)
- `core.py`: Changed `generate_report()` to write to `reports_dir` instead of `data_dir`

**Status:** RESOLVED

---

## Issue 8: `run_enhanced_analysis.py` modifies wrong config keys

**Diagnosis:** Confirmed. Script uses non-existent config keys.

**Action:** This script was already moved to `old/` in a prior session. It is no longer part of the active codebase.

**Status:** RESOLVED (previously archived)

---

## Issue 9: Effect-estimation layer is proxy logic, not real recovery

**Diagnosis:** Confirmed. The estimators assume wrong functional forms:
- `estimate_beta()` assumes linear ecDNA→VEGF, but simulator uses EGFR with sqrt scaling
- `estimate_delta()` assumes linear ecDNA→migration, but simulator uses EGFR with linear+hypoxia
- `estimate_alpha()` is a rough proxy from division event counts

**Action:**
- `estimate_beta()`: Now uses `sqrt(egfr_expression)` as the regressor, matching the simulator's `S_base * (1 + β * sqrt(EGFR))` form
- `estimate_delta()`: Now uses `egfr_expression` as the regressor and filters to normoxic cells to remove hypoxia confounding, matching `v_base * (1 + δ * EGFR)`
- `estimate_alpha()`: Now uses inter-division intervals from lineage data grouped by EGFR expression, fitting the `T_base / (1 + α * log2(1 + EGFR))` model
- All estimators now use EGFR as the exposure variable (consistent with Issue 1 fix)

**Status:** RESOLVED

---

## Issue 10: Tumor seeding doesn't clear existing cells

**Diagnosis:** Confirmed. Comment says "Clear any existing cell" but code doesn't.

**Action:**
- `initialization.py`: `seed_tumor()` now checks for and removes any existing cell at the seed position before placing the tumor cell

**Status:** RESOLVED

---

## Issue 11: Documentation has materially drifted from implementation

**Diagnosis:** Confirmed. Multiple contradictions across README, tutorial, analysis guide, and immune system docs.

**Action:**
- All documentation updated after code fixes to reflect the consistent causal chain: ecDNA → EGFR → phenotype
- Removed "reduced form" language from analysis_guide.md
- Updated immune_system.md parameter tables to match actual default.json values
- README causal DAG now shows EGFR as explicit exposure node

**Status:** RESOLVED

---

## Issue 12: Engineering bugs (CausalEdge templates, pandas dependency)

**Diagnosis:** Confirmed.
- `CausalEdge.equation` format only provides `{param}`, `{source}`, `{target}` but some templates use `{base}`
- Scripts import pandas but it's not in pyproject.toml dependencies

**Action:**
- `causal_dag.py`: Added `base="base"` to the format call as a safe default; templates that use `{base}` now render correctly
- `pyproject.toml`: Scripts using pandas were already moved to `old/`. No pandas dependency needed in active code.

**Status:** RESOLVED

---

## Issue 13: No automated test suite

**Diagnosis:** Confirmed. No tests/ directory exists.

**Action:**
- Created `tests/` directory with:
  - `test_segregation.py`: Tests ecDNA conservation, expected daughter fractions under replication model
  - `test_iv_pipeline.py`: Tests that D ≠ Z, EGFR column round-trips through loader
  - `test_discovery.py`: Tests that PC and GES produce acyclic graphs
  - `test_causal_chain.py`: Tests that all phenotype modulation routes through EGFR
  - `test_output_paths.py`: Tests that viewer HTML references correct asset paths

**Status:** RESOLVED

---

## Biological Best-Practice Gaps

**Diagnosis:** The review correctly identifies that CAUSANTA is a stylized benchmark, not a validated GBM simulator. The gaps (single tumor state, single ecDNA axis, light tissue physics, no calibration/UQ) are real but are design scope decisions, not bugs.

**Action:** Added a "Scope and Limitations" section to README.md and instructions.md explicitly positioning the project as a **causal inference benchmark prototype**, not a biologically faithful GBM simulator. The roadmap for biological extensions (tumor state programs, niche logic, calibration) is documented in instructions.md Priority 4-5.

**Status:** RESOLVED (documented as future work)

---

## Summary

| Issue | Severity | Status |
|-------|----------|--------|
| 1. IV benchmark D=Z | Critical | RESOLVED |
| 2. Exclusion restriction violated | Critical | RESOLVED |
| 3. Segregation math wrong | Critical | RESOLVED |
| 4. Unit-inconsistent physics | Major | RESOLVED (documented) |
| 5. GES produces cycles | Major | RESOLVED |
| 6. Immune chemotaxis dead | Major | RESOLVED |
| 7. Output path mismatch | Major | RESOLVED |
| 8. Wrong config keys in script | Major | RESOLVED (archived) |
| 9. Effect estimation proxies | Major | RESOLVED |
| 10. Tumor seeding no-clear | Minor | RESOLVED |
| 11. Documentation drift | Major | RESOLVED |
| 12. Engineering bugs | Minor | RESOLVED |
| 13. No test suite | Major | RESOLVED |
| Bio best-practice gaps | Scope | RESOLVED (documented) |
