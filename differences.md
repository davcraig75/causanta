# Differences Between Causal.v042826.v6.docx and Current Manuscript/Simulation

This document identifies discrepancies between Causal.v042826.v6.docx and the current state of the manuscript, code, and simulation results. **Items marked with [HIGHLIGHT] should be highlighted yellow** as they represent corrections needed from v6.

---

## 1. Structural Equation Parameters

### Gene Dosage (Z → X)

**v6.docx (INCORRECT):**
- X_base = 1.0
- κ = 0.5
- Formula: `X := X_base + κ · min(Z, Z_max) + ε_X`

**[HIGHLIGHT] Current/Correct:**
- **X_base = 2.89** (reflects population-level effects including hypoxia)
- **κ = 1.21** (expression contribution per ecDNA copy)
- Formula: `X := 2.89 + 1.21 · min(Z, Z_max) + ε_X`

*Rationale:* Parameters calibrated to match first-stage regression results from simulation. The normoxic cell regression confirms slope of 1.21.

---

### Hypoxia Threshold

**v6.docx (INCORRECT):**
- τ_hypoxia = 10 mmHg

**[HIGHLIGHT] Current/Correct:**
- **τ_hypoxia = 18 mmHg**

*Rationale:* Configured in paper_config_6mm.json as `hypoxia_threshold_mmHg: 18`.

---

### Cell Cycle Duration / Division Time

**v6.docx (INCORRECT):**
- T_base = 36 hours

**[HIGHLIGHT] Current/Correct:**
- **T_base = 24 hours**

*Rationale:* Configured in paper_config_6mm.json as `division_time_mean_hr: 24`.

---

### Apoptosis Effect (γ)

**v6.docx (INCORRECT):**
- γ = 0.2

**[HIGHLIGHT] Current/Correct:**
- **γ = 0.5**

*Rationale:* Configured in paper_config_6mm.json as `ecDNA_effect_on_survival: 0.5`.

---

### Oxygen Supply Parameters

**v6.docx (INCORRECT):**
- q = 5 h⁻¹
- O2_blood = 60 mmHg

**[HIGHLIGHT] Current/Correct:**
- **q = 2 h⁻¹**
- **O2_blood = 40 mmHg**

*Rationale:* Configured in paper_config_6mm.json.

---

### VEGF Secretion Equation

**v6.docx (INCORRECT):**
- Uses linear form: `(1 + β · X)`

**[HIGHLIGHT] Current/Correct:**
- Uses square-root form: `(1 + β · √X)`

*Rationale:* Square-root form reflects receptor saturation kinetics. See ecdna.py `modulate_vegf_secretion()`.

---

## 2. Table 2: Cell Type Parameters

**v6.docx (INCORRECT):**
- Tumor division time: 36h

**[HIGHLIGHT] Current/Correct:**
- **Tumor division time: 24h**

---

## 3. Simulation Statistics

### Division Events

**v6.docx:**
- 51,568 division events from 168-hour simulation

**[HIGHLIGHT] Current Simulation:**
- **51,658 division events from 300-hour simulation on 6mm × 6mm domain**

---

### Tumor Cell Count at Final Timepoint

**v6.docx:**
- 4,823 tumor cells

**[HIGHLIGHT] Current Simulation:**
- **44,400 tumor cells at t = 300h**

---

### First-Stage Regression

**v6.docx (INCORRECT):**
- EGFR = 1.02 + 0.49 · ecDNA
- R² = 0.87
- F-statistic = 3,420

**[HIGHLIGHT] Current/Correct:**
- **EGFR = 2.89 + 1.21 · ecDNA** (for normoxic cells; hypoxic cells show 2.5× HIF-1α upregulation)
- **R² = 0.851**
- **F-statistic = 253,699** (exceeds threshold by >25,000-fold)

---

### Text Describing κ

**v6.docx (INCORRECT):**
- "each additional ecDNA copy increases expression by approximately κ=0.5 units"

**[HIGHLIGHT] Current/Correct:**
- **"each additional ecDNA copy increases expression by approximately κ=1.21 units"**

---

## 4. OLS Bias Claims

**v6.docx:**
- "OLS overestimates the true causal effect by 140%"
- Reports 4 parameters (α, β, δ, γ) with specific estimates

**[HIGHLIGHT] Current/Correct:**
- **"OLS overestimates the true causal effect on VEGF secretion by 117% and on migration by 27%"**
- Focus on 2 parameters (β for VEGF, δ for migration) with updated estimates

---

## 5. Figure 3 / DAG Diagram

**v6.docx (INCORRECT):**
- Shows κ=0.5 in diagram
- Shows γ=0.2 in diagram

**[HIGHLIGHT] Current/Correct:**
- **κ=1.21** in diagram
- **γ=0.5** in diagram

---

## 6. Figure 1 Legend

**v6.docx (INCORRECT):**
- "51,658 division events from a 300-hour CAUSANTA simulation" (correct number, but variance = 0.0124)

**[HIGHLIGHT] Current/Correct:**
- Variance = **0.0071** (matches theoretical for Binomial(N, 0.5))

---

## 7. Minor Text Corrections

### Section: "Two-Stage Least Squares Estimation"

**v6.docx (INCORRECT):**
- Contains "wheretThe" (typo)

**[HIGHLIGHT] Current/Correct:**
- **"where the"** (corrected)

---

## Summary of All Yellow Highlights Needed

1. **Gene dosage equation:** Change X_base=1.0→2.89, κ=0.5→1.21
2. **Hypoxia threshold:** Change τ_hypoxia=10→18 mmHg
3. **Division time:** Change T_base=36→24 hours
4. **Apoptosis effect:** Change γ=0.2→0.5
5. **Oxygen supply:** Change q=5→2 h⁻¹, O2_blood=60→40 mmHg
6. **VEGF equation:** Change linear to sqrt form
7. **Table 2 Tumor division:** Change 36h→24h
8. **Division events:** Update to 51,658 from 300h simulation
9. **Tumor cells:** Update to 44,400
10. **First-stage regression:** Update to EGFR=2.89+1.21·ecDNA, R²=0.851, F=253,699
11. **κ in text:** Multiple locations, change 0.5→1.21
12. **DAG diagram:** Update κ=1.21, γ=0.5
13. **OLS bias claims:** Update percentages
14. **Segregation variance:** Update to 0.0071
15. **"wheretThe" typo:** Fix to "where the"

---

*Document generated by comparing Causal.v042826.v6.docx against current manuscript.md, ecdna.py, and paper_config_6mm.json.*
