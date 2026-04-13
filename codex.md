# CAUSANTA Deep Review

Date: 2026-04-13

## Scope

This report is based on a deep read of the repository's documentation and source code, plus light runtime validation.

Reviewed directly:

- Top-level docs and narrative: `README.md`, `Instructions.md`, `docs/tutorial.md`, `docs/analysis_guide.md`, `docs/immune_system.md`
- Simulation engine: everything under `causanta/simulate/`
- Analysis stack: everything under `causanta/analyze/`
- Ground-truth graph layer: `causanta/graph/causal_dag.py`
- User-facing scripts and examples: `scripts/`, `examples/`
- Output artifacts from existing runs under `output/`

Runtime checks performed:

- `python -m causanta.simulate.core --hours 5`
- `python examples/analyze_simulation.py output/run_20260413_093442/data`
- Direct probe of `discover_causal_structure(..., algorithm='ges')` on `output/run_20260413_093442/data`

External best-practice cross-checks used for context:

- PhysiCell as a reference architecture for multicellular cancer simulation: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005991
- Recent ABM calibration and uncertainty-quantification practice: https://pmc.ncbi.nlm.nih.gov/articles/PMC10869399/
- Recent GBM spatial niche evidence: https://link.springer.com/article/10.1186/s40478-024-01769-0
- Recent ecDNA spatial heterogeneity/evolution evidence in GBM: https://pmc.ncbi.nlm.nih.gov/articles/PMC12498097/

## Executive Verdict

CAUSANTA is a thoughtful and ambitious prototype for a synthetic spatial-causal benchmark. Its strongest ideas are clear:

- make the causal structure explicit
- generate tissue-like spatial data rather than tabular toy data
- use ecDNA segregation as a biologically motivated source of quasi-random variation
- pair simulation, reporting, and analysis in one repo

However, in its current form, it is **not yet using best practice for a rigorous glioma section simulation benchmark**, and it is **not yet a biologically trustworthy glioma simulator**.

The main reason is not lack of ambition. It is that the repo's three core contracts are currently out of sync:

1. the docs say the benchmark is `ecDNA -> EGFR -> phenotype`
2. parts of the simulator actually use `ecDNA -> phenotype` directly
3. the analyzer often treats `ecDNA` as both instrument and treatment

That breaks the central scientific promise of the repository: a known, internally consistent causal data-generating process that can be used to benchmark causal inference methods.

My bottom-line assessment:

| Area | Assessment |
|---|---|
| Overall code modularity | Good prototype quality |
| Causal benchmark validity | Currently compromised |
| Biological fidelity to GBM tissue sections | Low to moderate |
| Numerical architecture | Reasonable prototype, not calibrated |
| Analysis-method correctness | Several important flaws |
| Documentation consistency | Drifted and contradictory |
| Readiness for serious method benchmarking | Not yet |

The repo is closest to a **research prototype / educational benchmark**. It is not yet a robust foundation for publishable causal benchmarking or for claiming best-practice glioma section simulation.

## What This Repository Is Trying To Do

At a high level, CAUSANTA combines four systems:

1. A synthetic glioma tissue-section simulator.
2. A formal causal-DAG layer describing the intended ground truth.
3. An analysis layer that tries to recover effects and structure from generated data.
4. A reporting/visualization layer that packages the outputs.

The intended research story is:

- ecDNA copy number varies stochastically through mitosis
- ecDNA changes EGFR expression
- EGFR changes tumor behaviors such as proliferation, migration, VEGF secretion, and survival
- local oxygen and nutrient fields confound phenotype-outcome relationships
- causal methods should recover the true effects better than naive association methods

That is a strong concept. It is exactly the kind of benchmark design that could be valuable if the internal contract is made consistent and validated.

## Repository Structure And Architecture

### 1. Simulation Layer

The simulation lives under `causanta/simulate/`.

Core architectural pieces:

- `config.py`: frozen dataclass config loader from JSON
- `domain.py`: spatial coordinate transforms and interpolation
- `environment.py`: grid-based reaction-diffusion fields
- `cells.py`: cell state and spatial indexing
- `behaviors.py`: cell cycle, death, migration, immune interaction
- `ecdna.py`: ecDNA inheritance and phenotype modulation
- `initialization.py`: normal tissue, vasculature, tumor seeding
- `angiogenesis.py`: VEGF-driven vessel sprouting
- `io.py`: TSV/log output
- `viewer.py`, `visualization.py`, `reporting.py`: visualization/report generation

Simulation loop in `causanta/simulate/core.py`:

1. diffuse/update environment
2. process cell behaviors in shuffled order
3. run angiogenesis
4. recruit immune cells
5. clear necrotic cells and update occupancy
6. write outputs

This is a clean and understandable structure for a prototype. The modularity is a genuine strength.

### 2. Cell/Environment Representation

The code uses:

- continuous cell positions in a 2D plane
- a coarser regular environment grid
- per-cell state variables for lineage, ecDNA, hypoxia, migration, VEGF, immune state
- shape-aware collision checks for placement and movement

This is a reasonable architecture for a tissue-section simulator. It is more meaningful than a purely lattice-bound toy model and simpler than a full mechanical cell-body simulator.

### 3. Analysis Layer

The analysis lives under `causanta/analyze/`.

It includes:

- TSV loading utilities
- effect estimation
- IV / 2SLS routines
- causal discovery
- matching / IPW / regression
- SEM and sensitivity analysis
- HTML reporting

Architecturally, this is good in one sense: the benchmark includes the "recover the truth" side of the workflow. But scientifically, this layer is where some of the deepest problems currently appear.

### 4. Ground-Truth Graph Layer

`causanta/graph/causal_dag.py` formalizes the intended causal model. This is conceptually excellent. A simulation benchmark should have a first-class representation of its data-generating DAG.

The problem is that the DAG, the simulator, and the analyzer do not currently agree with each other.

## What The Repo Gets Right

These are real strengths, not just polite positives:

- The repo has a clear scientific identity rather than being a generic tumor toy.
- The simulation, analysis, and graph layers are separated sensibly.
- Configuration is explicit and reproducible via JSON plus RNG seed.
- The simulation writes useful benchmark outputs: cell snapshots, environment grids, lineage, summary logs.
- The code is readable. Important logic is not buried behind opaque abstractions.
- The repo tries to expose both biological assumptions and statistical assumptions.
- Shape-aware collision and continuous-space placement are nice touches for a 2D tissue-section engine.
- The reaction-diffusion solver is structured carefully enough to be numerically useful as a prototype.

If the core causal contract were repaired, this could become a genuinely useful benchmark platform.

## Best-Practice Assessment

### For A Causal-Inference Benchmark

Best practice for this repo's stated purpose would mean:

- the data-generating SCM is explicit
- the simulator actually implements that SCM
- the emitted data contain the variables needed to test the intended estimators
- the recovery methods are statistically valid for the generated data
- benchmark claims are stress-tested across seeds, parameter regimes, and assumption violations

CAUSANTA only partially satisfies this today.

It has:

- an explicit intended DAG
- reproducible generation
- relevant outputs

But it does not yet have:

- a simulator that consistently implements the stated DAG
- an analyzer that uses the stated exposure/instrument correctly
- a reliable mapping from configured "truth" to estimated effects

So for benchmark purposes, the current implementation is **conceptually promising but scientifically unstable**.

### For A Biologically Faithful Glioma Section Simulator

Best practice for a modern GBM section simulator would usually include most of the following:

- explicit calibration and uncertainty quantification
- niche-aware microenvironmental structure
- realistic state heterogeneity, not just one generic tumor cell state
- mechanistically consistent tissue physics or clearly declared reduced-form physics
- immune behavior tied to actual gradients/contacts rather than mostly descriptive parameters
- validation against known GBM spatial patterns such as perivascular, perinecrotic, and immunosuppressive niches

Relative to that standard, CAUSANTA is currently a **stylized synthetic model**, not a validated GBM simulator.

That is not inherently bad. It only becomes bad when the repo presents the model as more realistic or more statistically rigorous than it currently is.

## Critical Findings

## 1. The IV benchmark is currently invalid in the analyzer

This is the single most important problem in the repo.

The intended story is:

- `ecDNA_count` is the instrument
- `EGFR_expression` is the exposure
- migration / VEGF / proliferation are outcomes

The output writer actually stores `egfr_expression` in cell TSVs:

- `causanta/simulate/io.py:24-32`
- `causanta/simulate/io.py:75-100`

But the analysis loader drops it:

- `causanta/analyze/loader.py:67-95`

Then the IV estimator uses `ecDNA_count` as both instrument and treatment:

- `causanta/analyze/iv.py:191-214`

Specifically:

- `Z = ecDNA_count`
- `D = Z`

That means the "first stage" is effectively regressing a variable on itself.

Observed consequence in runtime:

- `python examples/analyze_simulation.py output/run_20260413_093442/data`
- reported first-stage F-statistic: `122777165769784676409409536.0`

That is not evidence of a spectacular instrument. It is a symptom of a broken setup.

Why this is fundamental:

- the benchmark can no longer claim to test whether IV recovers the true causal effect of EGFR
- the 2SLS estimate collapses toward a dressed-up reduced-form regression
- the repo's central README/tutorial claim becomes false in practice

Impact:

- invalidates the strongest causal benchmarking claim in the repo
- makes any "IV beats OLS" demo scientifically unreliable

## 2. The simulator violates its own stated exclusion-restriction story

The docs repeatedly say ecDNA affects phenotype only through gene expression:

- `README.md:103-109`
- `docs/tutorial.md:58-75`
- `docs/analysis_guide.md:62-78`

The code does not implement that consistently.

Migration and VEGF are modulated through EGFR:

- `causanta/simulate/behaviors.py:92-130`

But apoptosis still uses `cell.ecDNA_count` directly:

- `causanta/simulate/behaviors.py:168-175`

And cell-cycle entry from G0 uses `cell.ecDNA_count` directly:

- `causanta/simulate/behaviors.py:417-433`

Those calls are inconsistent with the intended function signatures in `ecdna.py`, which are written in terms of EGFR expression:

- `causanta/simulate/ecdna.py:155-187`
- `causanta/simulate/ecdna.py:273-304`

So the implemented simulator is a hybrid:

- `ecDNA -> EGFR -> migration / VEGF`
- `ecDNA -> division`
- `ecDNA -> survival`

That breaks the clean causal story and weakens the claimed IV logic.

It gets worse:

after division, the parent's next cycle time is recomputed from EGFR,
but the daughter's first cycle time is not:

- parent reset and recompute: `causanta/simulate/behaviors.py:520-532`
- daughter created in G1 with default sampled cycle time: `causanta/simulate/behaviors.py:497-513`

This means two daughters do not receive symmetric ecDNA-dependent proliferation dynamics.

Why this matters:

- the true data-generating process is not the DAG being advertised
- the exclusion restriction is no longer cleanly implemented
- the benchmark truth is not well defined

## 3. The segregation analysis is mathematically inconsistent with the implemented inheritance model

The inheritance code does:

1. replication during S phase
2. random segregation of the replicated copies

See:

- `causanta/simulate/ecdna.py:85-130`

This means if a parent starts a cycle with `N` ecDNA copies, each daughter should inherit roughly `N`, not `N/2`, when compared to the pre-replication parent count.

But the analysis code assumes the expected daughter fraction relative to `parent_ecDNA_before` is `0.5`:

- `causanta/analyze/effects.py:382-430`

The tutorial/docs make the same assumption:

- `README.md:105-109`
- `docs/tutorial.md:68-76`
- `docs/analysis_guide.md:92-99`

Runtime evidence from actual lineage output:

- `output/run_20260413_093442/data/lineage.tsv` begins with `20 -> 19 + 20`
- daughter fraction there is `20 / 20 = 1.0`, not `0.5`

So today:

- the simulation is modeling "replicate then segregate"
- the analysis and docs are testing "segregate pre-replication copies"

This is a deep contract error, not a cosmetic doc issue.

Impact:

- current segregation statistics are misinterpreted
- current "instrument validation" diagnostics are not aligned with the simulator
- tutorial explanations of ecDNA inheritance are wrong for this codebase

## 4. The microenvironment physics are unit-inconsistent and therefore not quantitatively trustworthy

The environment fields are labeled with physical units:

- O2 in mmHg
- glucose in mM
- VEGF in nM
- lactate in mM

See docs/config:

- `causanta/simulate/params/default.json`
- `causanta/simulate/environment.py:30-54`

But the reaction terms mix incompatible units directly. Example:

- vascular O2 source is computed as a transfer coefficient times `(O2_blood - env.O2)` in mmHg
- cellular O2 sink is computed from `O2_consumption_amol_hr / voxel_area`

See:

- `causanta/simulate/environment.py:292-311`

The same problem exists for glucose.

There is no conversion from attomoles per hour into mmHg or mM. So the field equations are not physically dimensionally consistent. They are effectively phenomenological knobs dressed in physical units.

That may be acceptable in a toy benchmark if stated plainly. It is not best practice if the model is presented as a realistic glioma microenvironment simulator.

Impact:

- parameter values cannot be interpreted literally
- calibration against experimental measurements is not currently well posed
- downstream claims about oxygen or glucose realism should be treated cautiously

## 5. The causal discovery module can produce cyclic graphs

The GES implementation greedily adds edges without any acyclicity check:

- `causanta/analyze/discovery.py:215-328`

Nothing in the forward phase stops it from creating `A -> B` and `B -> A`.

I confirmed this by running the discovery code on repo-generated output:

- discovered edges included `glucose_local -> O2_local` and `O2_local -> glucose_local`
- discovered edges also included `migration_rate -> VEGF_secretion` and `VEGF_secretion -> migration_rate`

That is not a DAG.

The default discovery variables are also only a partial, renamed subset of the intended graph:

- `causanta/analyze/discovery.py:331-355`

They do not include the explicit exposure variable `EGFR_expression`, and they do not line up cleanly with node names in `causanta/graph/causal_dag.py`.

Impact:

- discovery results are not valid DAGs
- benchmark comparisons to "ground truth" are not on a stable variable mapping
- structural metrics can be misleading

## Major Issues And Needed Improvements

## 6. The immune model is documented as richer than it is actually wired

The docs describe chemokine-like recruitment and chemotaxis in detail:

- `docs/immune_system.md:66-140`

The config also exposes chemokine parameters:

- `causanta/simulate/params/default.json`
- `causanta/simulate/config.py`

But the actual environment has no chemokine field:

- `causanta/simulate/environment.py:37-46`

There is a `compute_immune_chemotaxis()` function:

- `causanta/simulate/behaviors.py:321-359`

But the main migration function never calls it:

- `causanta/simulate/behaviors.py:547-625`

So immune-cell `chemotaxis_chemokine` is effectively dead in the current movement code.

Also, `max_kills_before_exhaustion` is described in docs and config, but not actually used in `check_immune_kill()`:

- comment claims at `causanta/simulate/behaviors.py:195-196`
- implementation only increments `exhaustion_level` by `exhaustion_per_kill`: `causanta/simulate/behaviors.py:251-256`

This means the immune module is best described as a simplified contact/activation prototype, not a biologically realistic GBM immune simulation.

## 7. Reporting and viewer output paths are inconsistent

At simulation end:

- `visualization.vg.json` is written to `data/`
- `index.html` is written to `reports/`

See:

- `causanta/simulate/core.py:197-245`
- `causanta/simulate/visualization.py:201-243`

But `index.html` looks for `visualization.vg.json` in its own directory:

- `causanta/simulate/visualization.py:233-236`

Observed output from a runtime check:

- `output/run_20260413_095023/data/visualization.vg.json`
- `output/run_20260413_095023/reports/index.html`

So the viewer path is broken.

A related mismatch:

- the generated report lands in `data/report.html`
- the README says reports should live under `reports/`

That is a real UX issue for anyone relying on the documented output structure.

## 8. `run_enhanced_analysis.py` modifies the wrong config keys

This script appears intended to tune long runs, but it sets keys the config system does not use:

- wrong time key: `scripts/run_enhanced_analysis.py:103-105`
- wrong tumor-seed keys: `scripts/run_enhanced_analysis.py:108-111`

The actual config schema expects:

- `output_interval_hr`, not `output_interval_hours`: `causanta/simulate/config.py:19-25`
- tumor seeds under `tumor_seeds`, not `initialization.tumor_seed_count` or `initialization.tumor_seed_ecDNA`: `causanta/simulate/config.py:140-152`

So the script's printed settings do not fully match what the simulation actually receives.

## 9. The effect-estimation layer is mostly proxy logic, not real recovery of configured truth

Examples:

- `estimate_alpha()` uses mean log ecDNA at division divided by base cycle time, explicitly as a rough proxy: `causanta/analyze/effects.py:62-135`
- `estimate_beta()` assumes a linear ecDNA relationship, but the simulator uses EGFR and sqrt scaling: `causanta/analyze/effects.py:139-219`, versus `causanta/simulate/ecdna.py:190-222`
- `estimate_delta()` assumes linear ecDNA effects on migration, but the simulator uses EGFR plus hypoxia multiplier: `causanta/analyze/effects.py:222-285`, versus `causanta/simulate/ecdna.py:225-270`
- `estimate_gamma_from_lineage()` infers survival from enrichment, but the simulator does not emit enough death-history detail to make this a real survival analysis: `causanta/analyze/effects.py:288-357`

The repo currently reports configured-vs-estimated values as if they are directly comparable. In many cases, they are not.

That makes the report layer more didactic than rigorous.

## 10. Tumor seeding and docs imply behaviors that are not actually enforced

`seed_tumor()` says it will clear any existing cell at the tumor position:

- `causanta/simulate/initialization.py:294-311`

But it does not clear anything. It just adds the tumor cell.

In practice this may or may not collide with an existing initialized cell, depending on random placement. The comment is simply not implemented.

This is not as severe as the IV issues, but it is a good example of contract drift between intent and code.

## 11. Documentation has materially drifted from the implementation

There are several contradictory narratives across the docs:

- `README.md` and `docs/tutorial.md` say the exposure is explicit EGFR and exclusion is clean.
- `docs/analysis_guide.md:144-162` then says CAUSANTA can use ecDNA as both instrument and exposure in reduced form.
- `docs/immune_system.md` correctly describes an immunosuppressive default at lines `29-37`, but its parameter tables later list much less suppressive defaults at `83-122` than the actual `default.json`.

This makes it hard to know which scientific contract a user is supposed to trust:

- the README story
- the analysis-guide reduced-form story
- or the actual simulator

For a benchmark repo, that ambiguity is costly.

## 12. There are a few smaller but real engineering-contract bugs

Two worth calling out:

- `CausalEdge.equation` formats templates with `param`, `source`, and `target`, but some templates in the DAG builder also reference `{base}`. That can break serialization paths that call `edge.to_dict()`.
  - formatter: `causanta/graph/causal_dag.py:76-83`
  - example template with `{base}`: `causanta/graph/causal_dag.py:430-436`
- The package dependencies do not include `pandas`, but helper scripts import it.
  - dependencies: `pyproject.toml:10-19`
  - examples in scripts: `scripts/run_full_analysis.py:218-223`, `scripts/run_enhanced_analysis.py:156-159`

These are not as severe as the IV or inheritance issues, but they are the kind of drift that makes a scientific repo feel less reliable.

## 13. There is no automated test suite for the core scientific contracts

I did not find a tests directory or test files in the repo scan.

Given the issues above, the missing tests are significant. This project needs automated checks for at least:

- ecDNA conservation/inheritance expectations
- instrument/exposure/output schema consistency
- analysis invariants
- output-path correctness
- discovery acyclicity
- config-key coverage in helper scripts

For a scientific benchmark, "no tests" is not just an engineering weakness; it is a scientific reproducibility weakness.

## Biological Best-Practice Gap

Even after fixing the internal logic problems, the current simulator would still be a stylized GBM benchmark rather than a best-practice GBM section simulator.

Key gaps relative to current GBM spatial literature and mature multicellular simulation practice:

### A. Too little tumor-state heterogeneity

Recent GBM spatial work emphasizes distinct niches and segregated cell states, including OPC-like, NPC-like, AC-like, and MES-like programs, with strong perivascular and perinecrotic structure.

See:

- Liu et al. 2024: perinecrotic regions are more immunosuppressive and perivascular regions are more pro-inflammatory, and tumor cell states are spatially segregated. https://link.springer.com/article/10.1186/s40478-024-01769-0

CAUSANTA currently has one generic tumor type with ecDNA-driven modifiers. That is enough for a simplified benchmark, but not enough for a state-of-the-art GBM section model.

### B. Too little ecDNA evolutionary richness

Recent GBM ecDNA work highlights oncogene-specific spatial heterogeneity and evolutionary behavior, including distinct EGFR ecDNA dynamics and ecDNA-driven heteroplasmy/clonality.

See:

- Noorani et al. 2024/2025: ecDNA random segregation plus fitness differences create oncogene-specific spatial heterogeneity in GBM. https://pmc.ncbi.nlm.nih.gov/articles/PMC12498097/

CAUSANTA currently uses one EGFR-like ecDNA axis. That is acceptable for a first benchmark, but it is not best-practice ecDNA modeling.

### C. Tissue physics are lighter than mature multicellular frameworks

A mature reference point like PhysiCell emphasizes:

- lattice-free mechanics
- microenvironment-dependent phenotypes
- cell death, motility, and mechanical interactions as first-class modules
- explicit microenvironment fields and immune examples that tie movement to actual chemoattractant gradients

Reference:

- https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005991

CAUSANTA partially overlaps with this design philosophy, but its mechanics are much lighter:

- endpoint collision rather than richer force-based interactions
- 2D direct section simulation rather than a clearer 3D-to-2D sectioning story
- immune migration not actually wired to the documented chemotactic logic

### D. No serious calibration / validation / UQ workflow yet

Modern ABM practice increasingly treats calibration, uncertainty, and validation as core model-development steps rather than optional extras.

Reference:

- Wang et al. 2024 on calibrating tumor ABMs with approximate Bayesian computation and uncertainty estimation: https://pmc.ncbi.nlm.nih.gov/articles/PMC10869399/

CAUSANTA currently has:

- hand-set parameters
- didactic effect comparisons
- no calibration pipeline
- no uncertainty quantification around simulator parameters
- no validation targets against measured GBM spatial data

For a pedagogical prototype that is fine.
For "best practice", it is not.

## What I Would Call This Model Today

The most accurate description is:

**a modular, readable, causally motivated synthetic spatial benchmark prototype for glioma-like tissue, with meaningful promise but several unresolved scientific-contract bugs**

I would not currently describe it as:

- a validated glioma tissue simulator
- a rigorous IV benchmark
- a best-practice spatial GBM modeling platform

## Recommended Roadmap

## Priority 0: Decide The Scientific Contract

This must happen before further feature growth.

Choose one of these two paths:

### Path A: Explicit exposure benchmark

Keep the intended story:

- `ecDNA_count` = instrument
- `EGFR_expression` = exposure
- phenotypes depend on EGFR

Then:

- make all phenotype modulation use EGFR, not raw ecDNA
- emit and load EGFR everywhere
- make IV stage 1 use EGFR
- update discovery variables to include EGFR

### Path B: Reduced-form ecDNA benchmark

If the real goal is to benchmark methods on reduced-form ecDNA effects:

- rewrite the DAG
- stop claiming exclusion through EGFR
- treat ecDNA as treatment, not instrument-plus-exposure
- rewrite tutorials accordingly

Right now the repo is trying to occupy both positions at once.

## Priority 1: Fix The Current Scientific Bugs

In order:

1. Repair the IV pipeline so `D != Z`.
2. Make division and survival modulation consistent with the stated exposure.
3. Fix daughter-cell cycle initialization after division.
4. Align segregation diagnostics with the replication model.
5. Enforce acyclicity in discovery.

Until those are fixed, the benchmark truth is not stable enough.

## Priority 2: Add Scientific Regression Tests

Minimum high-value tests:

- ecDNA inheritance expectations under replication-plus-segregation
- EGFR column round-trip from simulation output into analysis loader
- IV sanity test on a small synthetic dataset where true effect is known
- discovery algorithm must not return cycles
- output viewers must resolve their assets
- helper scripts must only set schema-recognized config keys

## Priority 3: Separate Benchmark Fidelity From Biological Fidelity

Create two explicit tiers:

### Tier 1: Benchmark mode

- simple but internally consistent
- fast
- designed for many seeds and parameter sweeps
- explicit assumption checks and falsification settings

### Tier 2: Biological realism mode

- richer GBM state heterogeneity
- niche-aware immune and vascular logic
- calibrated units or clearly dimensionless reduced-form fields
- external validation targets

That separation would keep the project honest and much easier to maintain.

## Priority 4: If Biological Realism Matters, Expand The GBM State Model

Strong candidates:

- GBM state programs: OPC-like, NPC-like, AC-like, MES-like
- explicit perivascular and perinecrotic niche logic
- chemokine field instead of placeholder comments
- macrophage / TAM polarization instead of one recruited-immune bucket
- anisotropic migration substrates for white-matter-like tracks
- more explicit de novo ecDNA gain/loss or oncogene-specific ecDNA classes

## Priority 5: Add Calibration And Ensemble Benchmarking

For this repo's intended use, one-off simulations are not enough.

Add:

- multi-seed experiment sweeps
- parameter registries for benchmark scenarios
- calibration notebooks or scripts
- uncertainty summaries for key outputs
- benchmark datasets with known pass/fail expectations for each estimator

## Final Assessment

The architecture is promising. The codebase is readable. The scientific idea is strong.

But the current implementation has several deep inconsistencies that prevent me from calling it best practice for glioma section simulation or for causal benchmarking:

- the benchmark's instrument/exposure logic is broken in the analyzer
- the simulator does not consistently implement the DAG it advertises
- the segregation analysis is inconsistent with the implemented inheritance model
- the environment physics are not unit-consistent
- parts of the immune, reporting, and discovery systems are only partially wired or incorrect

If the project goal is:

- **education / prototype exploration**: it is already useful
- **serious causal-method benchmarking**: it needs core scientific repairs first
- **biologically faithful GBM section simulation**: it needs both those repairs and a second wave of biological calibration and niche/state modeling

## Source Links

- PhysiCell: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005991
- ABM calibration with ABC: https://pmc.ncbi.nlm.nih.gov/articles/PMC10869399/
- GBM spatial transcriptomic niches: https://link.springer.com/article/10.1186/s40478-024-01769-0
- ecDNA spatial heterogeneity in GBM: https://pmc.ncbi.nlm.nih.gov/articles/PMC12498097/
