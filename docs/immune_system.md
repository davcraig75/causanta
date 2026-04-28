# CAUSANTA Immune System Module

## Overview

The immune system module in CAUSANTA simulates **cytotoxic T lymphocyte (CTL) and microglia-mediated tumor killing** using biologically realistic immunological synapse mechanics. Unlike simple probabilistic death rates, this implementation models the time-dependent formation of immunological synapses, CTL exhaustion, and activation dynamics based on experimental literature.

## Glioblastoma and Immune Evasion

**Glioblastoma multiforme (GBM) is characterized by profound immunosuppression.** The default parameters in CAUSANTA reflect this biology:

### Why Tumors Typically "Win" in GBM

1. **Blood-Brain Barrier (BBB)**: The CNS is an immune-privileged site with limited T cell access. Even when the BBB is disrupted by tumor growth, infiltration remains restricted.

2. **Tumor-Associated Macrophages/Microglia (TAMs)**: Up to 30-50% of the tumor mass in GBM consists of TAMs. However, these cells are often "re-educated" by the tumor to adopt an immunosuppressive M2-like phenotype, promoting tumor growth rather than killing.

3. **Immunosuppressive Cytokines**: GBM cells secrete TGF-beta, IL-10, and PD-L1, which:
   - Inhibit T cell activation and proliferation
   - Induce T regulatory cells (Tregs)
   - Promote T cell exhaustion

4. **Metabolic Reprogramming**: The hypoxic, acidic tumor microenvironment (Warburg effect) impairs T cell function:
   - Low glucose limits T cell glycolysis
   - High lactate inhibits T cell proliferation
   - Hypoxia stabilizes HIF-1alpha in tumor cells, enhancing survival

5. **Reduced MHC Expression**: GBM cells often downregulate MHC class I, becoming "invisible" to CTLs.

### Default Parameters Reflect This Biology

The default immune parameters create an immunosuppressive environment where:
- Kill probability is low (25% per synapse vs 70% in standard CTLs)
- Exhaustion occurs rapidly (after ~3 kills)
- Activation requires high tumor proximity
- Deactivation occurs faster than activation

This results in tumors that grow despite immune pressure, matching clinical observations.

## Biological Background

### Immunological Synapse Formation

When a CTL encounters a target tumor cell, it does not kill instantly. Instead:

1. **Recognition Phase (~30 min)**: The T cell receptor (TCR) recognizes MHC-presented antigens
2. **Synapse Formation (1-2 hours)**: A stable immunological synapse forms with organized signaling domains
3. **Cytotoxic Phase**: Perforin/granzyme release causes target cell death

This process is modeled as a contact-duration-dependent kill probability.

**Literature References:**
- Huppa JB, Davis MM. *T-cell-antigen recognition and the immunological synapse.* Nat Rev Immunol. 2003;3(12):973-83
- Dustin ML. *The immunological synapse.* Arthritis Res Ther. 2008;10(Suppl 1):S11

### Serial Killing and Exhaustion

CTLs can kill multiple targets sequentially ("serial killing"), but their efficiency declines after ~10 kills:

- **Initial Phase**: High killing efficiency
- **Exhaustion**: Progressive loss of cytotoxic granules
- **Terminal Phase**: Loss of proliferative capacity

**Literature Reference:**
- Boissonnas A, et al. *In vivo imaging of cytotoxic T cell infiltration and elimination of a solid tumor.* J Exp Med. 2007;204(2):345-56

### Chemotaxis

Immune cells are attracted to tumors via chemokine gradients:
- **CCL2/MCP-1**: Monocyte chemoattractant protein
- **CCL5/RANTES**: Recruits T cells and macrophages
- **CXCL10/IP-10**: IFN-gamma-induced chemokine

In CAUSANTA, chemotaxis is approximated using VEGF gradients (correlating with tumor hypoxia/metabolism).

---

## Configuration Parameters

All parameters are set in the `immune_recruitment` section of the JSON configuration file.

### Recruitment Parameters

| Parameter | Default | Unit | Description |
|-----------|---------|------|-------------|
| `chemokine_threshold` | 2.0 | nM | Minimum chemokine (VEGF proxy) for immune recruitment |
| `recruitment_rate_per_hr` | 0.003 | 1/hr | Rate of new immune cell entry at vascular sites (GBM-restricted) |

**How to adjust:**
- Increase `recruitment_rate_per_hr` to simulate stronger immune response (e.g., immunotherapy)
- Decrease to simulate immunosuppressive tumor microenvironment

### Killing Mechanics

| Parameter | Default | Unit | Description |
|-----------|---------|------|-------------|
| `kill_radius_um` | 15.0 | um | Distance to detect potential targets |
| `synapse_formation_time_hr` | 2.5 | hr | Time for immunological synapse to form (GBM: slower) |
| `kill_probability_per_synapse` | 0.25 | probability | P(kill) once synapse is complete (GBM: immunosuppressed) |
| `min_contact_for_kill_hr` | 1.5 | hr | Minimum contact before any kill possible |

**The killing process:**
```
If contact_duration < min_contact_for_kill:
    kill_probability = 0
Else:
    synapse_progress = contact_duration / synapse_formation_time
    base_kill_prob = kill_probability_per_synapse * synapse_progress
    effective_prob = base_kill_prob * (1 - exhaustion_penalty) * activation_level
```

**How to adjust:**
- Decrease `synapse_formation_time_hr` to model more efficient CTLs (e.g., CAR-T cells)
- Increase `kill_probability_per_synapse` for more potent effector cells

### Exhaustion Mechanics

| Parameter | Default | Unit | Description |
|-----------|---------|------|-------------|
| `max_kills_before_exhaustion` | 3 | count | CTLs exhaust after ~3 serial kills (GBM: rapid exhaustion) |
| `exhaustion_per_kill` | 0.35 | fraction | Exhaustion increment per kill (0.35 = ~3 kills to full exhaustion) |
| `exhaustion_recovery_rate_per_hr` | 0.005 | 1/hr | Very slow recovery (GBM: sustained suppression) |
| `exhausted_kill_penalty` | 0.8 | fraction | Kill probability multiplier when fully exhausted |

**Exhaustion model:**
```python
# On successful kill:
exhaustion_level += exhaustion_per_kill
exhaustion_level = min(1.0, exhaustion_level)

# Kill probability reduction:
effective_prob = base_prob * (1.0 - exhaustion_level * exhausted_kill_penalty)

# Recovery when not in contact:
if not in_contact:
    exhaustion_level -= exhaustion_recovery_rate_per_hr * dt
```

**How to adjust:**
- Decrease `exhaustion_per_kill` to model more persistent CTLs
- Increase `exhaustion_recovery_rate_per_hr` to model faster recovery (e.g., with IL-2 support)

### Activation Dynamics

| Parameter | Default | Unit | Description |
|-----------|---------|------|-------------|
| `activation_radius_um` | 30.0 | um | Distance to sense tumor for activation |
| `activation_rate_per_hr` | 0.1 | 1/hr | Rate of activation when tumor nearby |
| `deactivation_rate_per_hr` | 0.15 | 1/hr | Rate of deactivation when away from tumor (GBM: faster than activation) |
| `min_activation_for_kill` | 0.5 | level | Minimum activation to attempt killing |

**Activation model:**
- Immune cells start with `activation_level = 0.5` when recruited
- Activation increases when tumor cells are within `activation_radius_um`
- Activation decays when no tumor nearby
- Killing requires `activation_level >= min_activation_for_kill`

**How to adjust:**
- Increase `activation_rate_per_hr` to model pre-activated cells (e.g., ex vivo expanded TILs)
- Decrease `min_activation_for_kill` for more responsive immune cells

### Chemokine Diffusion (Advanced)

| Parameter | Default | Unit | Description |
|-----------|---------|------|-------------|
| `chemokine_diffusion_um2_hr` | 1800.0 | um^2/hr | Chemokine diffusion coefficient |
| `chemokine_decay_per_hr` | 0.5 | 1/hr | Chemokine decay rate |
| `tumor_chemokine_secretion` | 0.3 | relative | Tumor secretion rate |

**Note:** Currently chemokine is approximated using VEGF. These parameters prepare for future explicit chemokine substrate.

---

## Cell Types Involved

### Recruited Immune Cells (type 7: RecruitedImmune)

Represents infiltrating CTLs, NK cells, and macrophages:
- Enter from vasculature in response to chemokine signals
- High migration speed (35 um/hr) with chemotaxis toward tumor
- Can kill tumor cells after synapse formation
- Subject to exhaustion after serial kills

### Microglia (type 3: Microglia)

Resident brain immune cells:
- Already present in tissue at initialization
- Require activation (via `is_reactive` flag) before killing
- Lower killing efficiency than recruited immune cells
- Activated by chemokine signals near tumor

---

## State Tracking

Each immune cell tracks:

| Field | Type | Description |
|-------|------|-------------|
| `contact_target_id` | int | ID of tumor cell in contact (-1 if none) |
| `contact_duration_hr` | float | Time in contact with current target |
| `kills_performed` | int | Total successful kills |
| `exhaustion_level` | float | 0.0 (fresh) to 1.0 (exhausted) |
| `activation_level` | float | 0.0 (inactive) to 1.0 (fully activated) |
| `is_reactive` | bool | Whether cell is in reactive/activated state |

---

## Simulation Workflow

Each simulation step:

1. **Update Activation**: Check for nearby tumor, increase/decrease activation
2. **Update Exhaustion**: If not in contact, slowly recover exhaustion
3. **Find Targets**: Search within `kill_radius_um` for tumor cells
4. **Attempt Kill**: 
   - Update contact tracking (same target = increment duration, new target = reset)
   - Check if contact duration exceeds minimum
   - Calculate kill probability based on synapse progress, exhaustion, activation
   - If kill succeeds: remove tumor cell, increment exhaustion, reset contact

---

## Example Scenarios

### Standard Immunotherapy
```json
"immune_recruitment": {
    "recruitment_rate_per_hr": 0.02,
    "kill_probability_per_synapse": 0.8,
    "exhaustion_per_kill": 0.08
}
```

### CAR-T Cell Therapy
```json
"immune_recruitment": {
    "synapse_formation_time_hr": 0.5,
    "kill_probability_per_synapse": 0.9,
    "max_kills_before_exhaustion": 20
}
```

### Immunosuppressive Tumor
```json
"immune_recruitment": {
    "recruitment_rate_per_hr": 0.002,
    "activation_rate_per_hr": 0.1,
    "exhaustion_per_kill": 0.15
}
```

### Checkpoint Inhibitor Response
```json
"immune_recruitment": {
    "exhaustion_recovery_rate_per_hr": 0.1,
    "exhausted_kill_penalty": 0.2,
    "activation_rate_per_hr": 0.5
}
```

---

## Output Analysis

The simulation outputs track immune cell behavior:

### Cell TSV Files
Each `cells_t*.tsv` includes immune cell fields:
- `activation_level`: Current activation state
- `exhaustion_level`: Current exhaustion state
- `kills_performed`: Total kills by this cell

### Summary Log
The `summary.log` includes:
- `immune_cells`: Count of recruited immune cells

### Analyzing Immune Response

```python
import pandas as pd
import matplotlib.pyplot as plt

# Load cell data
cells = pd.read_csv('cells_t000100.tsv', sep='\t')

# Filter immune cells
immune = cells[cells['cell_type_name'] == 'RecruitedImmune']

# Plot exhaustion distribution
plt.hist(immune['exhaustion_level'], bins=20)
plt.xlabel('Exhaustion Level')
plt.ylabel('Count')
plt.title('CTL Exhaustion Distribution at t=100h')
```

---

## Calibration Notes

### Matching Experimental Data

If you have experimental data on:

1. **Killing rates**: Adjust `kill_probability_per_synapse` and `synapse_formation_time_hr`
2. **Infiltration density**: Adjust `recruitment_rate_per_hr`
3. **Exhaustion markers (PD-1, TIM-3)**: Use `exhaustion_level` as proxy
4. **Treatment response**: Modify parameters to match observed tumor regression

### Sensitivity Analysis

Key parameters to vary in sensitivity analysis:
1. `synapse_formation_time_hr`: Controls killing speed
2. `exhaustion_per_kill`: Controls immune durability
3. `recruitment_rate_per_hr`: Controls immune pressure
4. `kill_probability_per_synapse`: Controls killing efficiency

---

## References

1. Huppa JB, Davis MM. T-cell-antigen recognition and the immunological synapse. Nat Rev Immunol. 2003;3(12):973-83.

2. Dustin ML. The immunological synapse. Arthritis Res Ther. 2008;10(Suppl 1):S11.

3. Boissonnas A, et al. In vivo imaging of cytotoxic T cell infiltration and elimination of a solid tumor. J Exp Med. 2007;204(2):345-56.

4. Weigelin B, et al. Intravital third harmonic generation microscopy of collective melanoma cell invasion. Intravital. 2012;1(1):32-43.

5. Breart B, et al. Two-photon imaging of intratumoral CD8+ T cell cytotoxic activity during adoptive T cell therapy in mice. J Clin Invest. 2008;118(4):1390-7.

6. Wiedemann A, et al. Cytotoxic T lymphocytes kill multiple targets simultaneously via spatiotemporal uncoupling of lytic and stimulatory synapses. PNAS. 2006;103(29):10985-90.
