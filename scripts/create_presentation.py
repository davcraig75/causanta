#!/usr/bin/env python3
"""Create CAUSANTA overview PowerPoint presentation."""

from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.enum.shapes import MSO_SHAPE
from pptx.oxml.ns import nsmap
from pathlib import Path

# Color scheme (Nature Methods inspired)
PRIMARY = RGBColor(0x00, 0x77, 0xBB)  # Blue
SECONDARY = RGBColor(0x33, 0x33, 0x33)  # Dark gray
ACCENT = RGBColor(0xEE, 0x77, 0x33)  # Orange
LIGHT_BG = RGBColor(0xF8, 0xF9, 0xFA)  # Light gray


def add_title_slide(prs, title, subtitle=""):
    """Add a title slide."""
    slide_layout = prs.slide_layouts[6]  # Blank
    slide = prs.slides.add_slide(slide_layout)

    # Title
    left = Inches(0.5)
    top = Inches(2.5)
    width = Inches(9)
    height = Inches(1.5)

    title_box = slide.shapes.add_textbox(left, top, width, height)
    tf = title_box.text_frame
    p = tf.paragraphs[0]
    p.text = title
    p.font.size = Pt(44)
    p.font.bold = True
    p.font.color.rgb = PRIMARY
    p.alignment = PP_ALIGN.CENTER

    if subtitle:
        top = Inches(4.2)
        height = Inches(1)
        sub_box = slide.shapes.add_textbox(left, top, width, height)
        tf = sub_box.text_frame
        p = tf.paragraphs[0]
        p.text = subtitle
        p.font.size = Pt(24)
        p.font.color.rgb = SECONDARY
        p.alignment = PP_ALIGN.CENTER

    return slide


def add_section_slide(prs, title):
    """Add a section divider slide."""
    slide_layout = prs.slide_layouts[6]  # Blank
    slide = prs.slides.add_slide(slide_layout)

    # Background shape
    left = Inches(0)
    top = Inches(2.8)
    width = Inches(10)
    height = Inches(1.5)

    shape = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, left, top, width, height)
    shape.fill.solid()
    shape.fill.fore_color.rgb = PRIMARY
    shape.line.fill.background()

    # Title
    title_box = slide.shapes.add_textbox(Inches(0.5), Inches(3), Inches(9), Inches(1))
    tf = title_box.text_frame
    p = tf.paragraphs[0]
    p.text = title
    p.font.size = Pt(40)
    p.font.bold = True
    p.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)
    p.alignment = PP_ALIGN.CENTER

    return slide


def add_content_slide(prs, title, bullets, notes=None):
    """Add a content slide with bullet points."""
    slide_layout = prs.slide_layouts[6]  # Blank
    slide = prs.slides.add_slide(slide_layout)

    # Title
    title_box = slide.shapes.add_textbox(Inches(0.5), Inches(0.3), Inches(9), Inches(0.8))
    tf = title_box.text_frame
    p = tf.paragraphs[0]
    p.text = title
    p.font.size = Pt(32)
    p.font.bold = True
    p.font.color.rgb = PRIMARY

    # Underline
    line = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(0.5), Inches(1.1), Inches(9), Inches(0.03))
    line.fill.solid()
    line.fill.fore_color.rgb = PRIMARY
    line.line.fill.background()

    # Bullets
    content_box = slide.shapes.add_textbox(Inches(0.5), Inches(1.4), Inches(9), Inches(5.5))
    tf = content_box.text_frame
    tf.word_wrap = True

    for i, bullet in enumerate(bullets):
        if i == 0:
            p = tf.paragraphs[0]
        else:
            p = tf.add_paragraph()

        # Handle nested bullets
        if isinstance(bullet, tuple):
            text, level = bullet
            p.text = text
            p.level = level
        else:
            p.text = bullet
            p.level = 0

        p.font.size = Pt(20) if p.level == 0 else Pt(18)
        p.font.color.rgb = SECONDARY
        p.space_before = Pt(8)
        p.space_after = Pt(4)

    if notes:
        notes_slide = slide.notes_slide
        notes_slide.notes_text_frame.text = notes

    return slide


def add_two_column_slide(prs, title, left_title, left_bullets, right_title, right_bullets):
    """Add a two-column slide."""
    slide_layout = prs.slide_layouts[6]  # Blank
    slide = prs.slides.add_slide(slide_layout)

    # Title
    title_box = slide.shapes.add_textbox(Inches(0.5), Inches(0.3), Inches(9), Inches(0.8))
    tf = title_box.text_frame
    p = tf.paragraphs[0]
    p.text = title
    p.font.size = Pt(32)
    p.font.bold = True
    p.font.color.rgb = PRIMARY

    # Left column title
    left_title_box = slide.shapes.add_textbox(Inches(0.5), Inches(1.2), Inches(4.2), Inches(0.5))
    tf = left_title_box.text_frame
    p = tf.paragraphs[0]
    p.text = left_title
    p.font.size = Pt(22)
    p.font.bold = True
    p.font.color.rgb = ACCENT

    # Left column content
    left_box = slide.shapes.add_textbox(Inches(0.5), Inches(1.7), Inches(4.2), Inches(5))
    tf = left_box.text_frame
    tf.word_wrap = True
    for i, bullet in enumerate(left_bullets):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.text = bullet
        p.font.size = Pt(16)
        p.font.color.rgb = SECONDARY
        p.space_before = Pt(6)

    # Right column title
    right_title_box = slide.shapes.add_textbox(Inches(5.3), Inches(1.2), Inches(4.2), Inches(0.5))
    tf = right_title_box.text_frame
    p = tf.paragraphs[0]
    p.text = right_title
    p.font.size = Pt(22)
    p.font.bold = True
    p.font.color.rgb = ACCENT

    # Right column content
    right_box = slide.shapes.add_textbox(Inches(5.3), Inches(1.7), Inches(4.2), Inches(5))
    tf = right_box.text_frame
    tf.word_wrap = True
    for i, bullet in enumerate(right_bullets):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.text = bullet
        p.font.size = Pt(16)
        p.font.color.rgb = SECONDARY
        p.space_before = Pt(6)

    return slide


def add_code_slide(prs, title, code_text, description=""):
    """Add a slide with code block."""
    slide_layout = prs.slide_layouts[6]  # Blank
    slide = prs.slides.add_slide(slide_layout)

    # Title
    title_box = slide.shapes.add_textbox(Inches(0.5), Inches(0.3), Inches(9), Inches(0.8))
    tf = title_box.text_frame
    p = tf.paragraphs[0]
    p.text = title
    p.font.size = Pt(32)
    p.font.bold = True
    p.font.color.rgb = PRIMARY

    if description:
        desc_box = slide.shapes.add_textbox(Inches(0.5), Inches(1.1), Inches(9), Inches(0.6))
        tf = desc_box.text_frame
        p = tf.paragraphs[0]
        p.text = description
        p.font.size = Pt(18)
        p.font.color.rgb = SECONDARY
        code_top = 1.8
    else:
        code_top = 1.3

    # Code background
    code_bg = slide.shapes.add_shape(
        MSO_SHAPE.ROUNDED_RECTANGLE,
        Inches(0.5), Inches(code_top), Inches(9), Inches(5.5 - (code_top - 1.3))
    )
    code_bg.fill.solid()
    code_bg.fill.fore_color.rgb = RGBColor(0x2D, 0x2D, 0x2D)
    code_bg.line.fill.background()

    # Code text
    code_box = slide.shapes.add_textbox(Inches(0.7), Inches(code_top + 0.2), Inches(8.6), Inches(5))
    tf = code_box.text_frame
    tf.word_wrap = True
    p = tf.paragraphs[0]
    p.text = code_text
    p.font.size = Pt(12)
    p.font.name = "Courier New"
    p.font.color.rgb = RGBColor(0xE0, 0xE0, 0xE0)

    return slide


def add_table_slide(prs, title, headers, rows):
    """Add a slide with a table."""
    slide_layout = prs.slide_layouts[6]  # Blank
    slide = prs.slides.add_slide(slide_layout)

    # Title
    title_box = slide.shapes.add_textbox(Inches(0.5), Inches(0.3), Inches(9), Inches(0.8))
    tf = title_box.text_frame
    p = tf.paragraphs[0]
    p.text = title
    p.font.size = Pt(32)
    p.font.bold = True
    p.font.color.rgb = PRIMARY

    # Table
    num_rows = len(rows) + 1
    num_cols = len(headers)

    table = slide.shapes.add_table(
        num_rows, num_cols,
        Inches(0.5), Inches(1.3),
        Inches(9), Inches(0.4 * num_rows)
    ).table

    # Set column widths
    col_width = Inches(9 / num_cols)
    for i in range(num_cols):
        table.columns[i].width = col_width

    # Header row
    for i, header in enumerate(headers):
        cell = table.cell(0, i)
        cell.text = header
        cell.fill.solid()
        cell.fill.fore_color.rgb = PRIMARY
        p = cell.text_frame.paragraphs[0]
        p.font.size = Pt(14)
        p.font.bold = True
        p.font.color.rgb = RGBColor(0xFF, 0xFF, 0xFF)
        p.alignment = PP_ALIGN.CENTER

    # Data rows
    for row_idx, row in enumerate(rows):
        for col_idx, value in enumerate(row):
            cell = table.cell(row_idx + 1, col_idx)
            cell.text = str(value)
            p = cell.text_frame.paragraphs[0]
            p.font.size = Pt(12)
            p.font.color.rgb = SECONDARY
            if row_idx % 2 == 1:
                cell.fill.solid()
                cell.fill.fore_color.rgb = LIGHT_BG

    return slide


def create_presentation():
    """Create the full CAUSANTA presentation."""
    prs = Presentation()
    prs.slide_width = Inches(10)
    prs.slide_height = Inches(7.5)

    # ==================== TITLE ====================
    add_title_slide(
        prs,
        "CAUSANTA",
        "Causal Analysis Using Somatic And Neighborhood Tissue Architecture"
    )

    # ==================== OVERVIEW ====================
    add_content_slide(prs, "Executive Summary", [
        "Simulation framework for benchmarking causal inference in spatial biology",
        "Generates synthetic tumor tissue with known ground-truth causal structure",
        "Key innovation: ecDNA segregation as natural randomization mechanism",
        ("Satisfies instrumental variable assumptions", 1),
        ("Enables causal effect identification from observational data", 1),
        "Validates whether causal discovery methods recover true effects",
    ])

    add_content_slide(prs, "The Problem: Causal Inference in Biology", [
        "Randomized experiments are often infeasible",
        ("Cannot randomly assign mutations to patients", 1),
        "Observational data has confounding",
        ("Hypoxia affects both gene expression AND cell behavior", 1),
        "Ground truth is rarely known",
        ("No way to validate causal effect estimates", 1),
        "CAUSANTA solves this with embedded ground truth + natural randomization",
    ])

    # ==================== SCIENTIFIC BACKGROUND ====================
    add_section_slide(prs, "Scientific Background")

    add_content_slide(prs, "Extrachromosomal DNA (ecDNA) Biology", [
        "Circular DNA molecules found in ~30% of cancers",
        "Lack centromeres - cannot attach to spindle fibers",
        "Segregate randomly during mitosis: Binomial(N, 0.5)",
        "Carry oncogenes (EGFR, MYC, etc.)",
        "Create tumor heterogeneity through random inheritance",
        "This randomness is the foundation for our IV approach",
    ])

    add_content_slide(prs, "The Confounding Problem", [
        "In tumor tissue, hypoxia (low O2) creates confounding:",
        ("Hypoxic cells upregulate EGFR expression", 1),
        ("Hypoxic cells also proliferate faster (HIF pathway)", 1),
        "Naive regression conflates:",
        ("Direct EGFR → Proliferation effect", 1),
        ("Confounded path through O2", 1),
        "Need instrumental variable to isolate causal effect",
    ])

    # ==================== SIV FRAMEWORK ====================
    add_section_slide(prs, "Somatic Instrumental Variable (SIV) Framework")

    add_content_slide(prs, "ecDNA as Instrumental Variable", [
        "ecDNA copy number (Z) instruments for EGFR expression (D)",
        "Three IV assumptions satisfied:",
        ("Relevance: ecDNA carries EGFR → gene dosage effect", 1),
        ("Independence: Binomial(N, 0.5) segregation is random", 1),
        ("Exclusion: ecDNA affects phenotypes via gene expression only", 1),
        "This enables Two-Stage Least Squares (2SLS) estimation",
    ])

    add_code_slide(prs, "Causal DAG", """
         ecDNA (Instrument Z)
              |
              v
         EGFR expression (Treatment D)  <--  O2 (Confounder U)
              |                               |
              v                               v
         Proliferation (Outcome Y)  <---------+


    Stage 1: EGFR_i = pi_0 + pi_1 * ecDNA_i + noise
    Stage 2: Prolif_i = beta_0 + beta_1 * EGFR_hat_i + error

    beta_1 = causal effect (purged of confounding)
""", "The causal structure embedded in CAUSANTA")

    add_table_slide(prs, "IV Assumptions & Satisfaction",
        ["Assumption", "Requirement", "How ecDNA Satisfies"],
        [
            ["Relevance", "Z affects D", "Gene dosage: more ecDNA → more EGFR"],
            ["Independence", "Z ⊥ U", "Binomial segregation is random"],
            ["Exclusion", "Z→Y only via D", "ecDNA affects phenotypes through expression"],
        ]
    )

    # ==================== ARCHITECTURE ====================
    add_section_slide(prs, "System Architecture")

    add_two_column_slide(prs, "Three Main Subsystems",
        "SIMULATE",
        [
            "• Agent-based tumor growth",
            "• ecDNA dynamics & segregation",
            "• Environment PDEs (O2, glucose, VEGF)",
            "• Immune system modeling",
            "• 6-phase simulation loop",
        ],
        "ANALYZE",
        [
            "• IV/2SLS estimation",
            "• Bootstrap confidence intervals",
            "• Power analysis",
            "• Heterogeneity by region/hypoxia",
            "• Sensitivity analysis (E-value, Rosenbaum)",
        ]
    )

    add_content_slide(prs, "Simulation Loop (6 Phases)", [
        "Phase 1: Environment Diffusion",
        ("Solve reaction-diffusion PDEs via implicit LOD", 1),
        "Phase 2: Cell Behaviors (random order)",
        ("Division with ecDNA segregation, migration, death", 1),
        "Phase 3: Angiogenesis",
        ("VEGF-driven vessel sprouting", 1),
        "Phase 4: Immune Recruitment",
        ("T-cell recruitment based on tumor burden", 1),
        "Phase 5: Cleanup - Remove dead cells",
        "Phase 6: Output - Write cell states, environment, lineage",
    ])

    # ==================== CELL MODEL ====================
    add_section_slide(prs, "Cell Agent Model")

    add_table_slide(prs, "Cell Types",
        ["Type", "Division", "Migration", "ecDNA", "Role"],
        [
            ["Tumor", "Yes (36h)", "10 μm/hr", "Yes", "Primary subject"],
            ["Neuron", "No", "0", "No", "Normal tissue"],
            ["Astrocyte", "Slow (100h)", "3 μm/hr", "No", "Glial support"],
            ["Microglia", "No", "30 μm/hr", "No", "Resident immune"],
            ["Endothelial", "Slow", "10 μm/hr", "No", "Blood vessels"],
            ["Recruited Immune", "No", "30 μm/hr", "No", "T-cells"],
            ["Necrotic", "No", "0", "No", "Dead cells"],
        ]
    )

    add_content_slide(prs, "Cell State Vector", [
        "Identity: cell_id, parent_id, cell_type",
        "Spatial: x, y (integer μm), angle, shape_path",
        "Proliferation: cell_cycle_phase (G0/G1/S/G2/M), cycle_clock_hr",
        "Genomic: ecDNA_count, ecDNA_cargo (gene list)",
        "Environment: O2_local, glucose_local, is_hypoxic",
        "Lineage: generation, time_born_hr",
        "Effective rates: migration_rate, VEGF_secretion",
    ])

    # ==================== ENVIRONMENT ====================
    add_section_slide(prs, "Environment Field Model")

    add_content_slide(prs, "Reaction-Diffusion PDEs", [
        "General form: ∂E/∂t = D∇²E - λE + sources - sinks",
        "Primary fields:",
        ("O2 (mmHg): D = 6×10⁶ μm²/hr, supplied by vessels", 1),
        ("Glucose (mM): D = 2.4×10⁶ μm²/hr, metabolic fuel", 1),
        ("VEGF (nM): D = 3.6×10⁴ μm²/hr, angiogenesis trigger", 1),
        ("Lactate (mM): Warburg effect, pH proxy", 1),
        "Solved via implicit LOD (Thomas algorithm)",
        "Unconditionally stable - no CFL constraint",
    ])

    add_table_slide(prs, "Environment Fields",
        ["Field", "Units", "Diffusion (μm²/hr)", "Normal Value", "Role"],
        [
            ["O2", "mmHg", "6.0 × 10⁶", "38", "Proliferation, necrosis"],
            ["Glucose", "mM", "2.4 × 10⁶", "5.0", "Metabolic fuel"],
            ["VEGF", "nM", "3.6 × 10⁴", "0", "Angiogenesis trigger"],
            ["Lactate", "mM", "9.0 × 10⁵", "1.0", "pH proxy"],
        ]
    )

    # ==================== ECDNA MECHANISM ====================
    add_section_slide(prs, "ecDNA Segregation (The SIV Mechanism)")

    add_code_slide(prs, "ecDNA Segregation During Division", """
    # During cell division (M phase completion):

    parent_ecDNA = N  # Current copy number

    # Binomial segregation - THE KEY MECHANISM
    daughter_ecDNA = Binomial(N, p=0.5)
    parent_ecDNA_after = N - daughter_ecDNA

    # This creates random variation INDEPENDENT of environment
    # Sisters in same microenvironment have different ecDNA
    # This mimics randomization within strata
""", "The core mechanism enabling causal inference")

    add_content_slide(prs, "ecDNA Modulation of Cell Behavior", [
        "Ground-truth causal effects embedded in simulation:",
        "Division: T_eff = T_base / (1 + α·log₂(1 + ecDNA))",
        ("α = 0.30: 30% faster division per ecDNA doubling", 1),
        "VEGF: S_eff = S_base × (1 + β·ecDNA)",
        ("β = 0.10: 10% more VEGF per ecDNA copy", 1),
        "Migration: v_eff = v_base × (1 + δ·ecDNA)",
        ("δ = 0.05: 5% faster migration per ecDNA copy", 1),
        "Survival: a_eff = a_base / (1 + γ·ecDNA)",
        ("γ = 0.20: 20% survival boost per ecDNA copy", 1),
    ])

    # ==================== ANALYSIS ====================
    add_section_slide(prs, "Causal Analysis Methods")

    add_content_slide(prs, "Two-Stage Least Squares (2SLS)", [
        "Stage 1 (First Stage): Predict treatment from instrument",
        ("EGFR_i = π₀ + π₁·ecDNA_i + η_i", 1),
        ("Get predicted values: EGFR_hat", 1),
        "Stage 2 (Second Stage): Regress outcome on prediction",
        ("Y_i = β₀ + β₁·EGFR_hat_i + ε_i", 1),
        "β₁ is the causal effect, purged of confounding",
        "Compare IV vs OLS to quantify confounding bias",
    ])

    add_table_slide(prs, "IV Diagnostics",
        ["Diagnostic", "Purpose", "Threshold"],
        [
            ["F-statistic", "Instrument strength", "> 10 (strong)"],
            ["Wu-Hausman", "Endogeneity detection", "p < 0.05 → use IV"],
            ["Anderson-Rubin", "Weak-IV robust CI", "Valid with weak instruments"],
            ["Rosenbaum bounds", "Sensitivity to confounding", "Report Γ*"],
            ["E-value", "Min confounding to explain", "> 2.0 (robust)"],
        ]
    )

    add_content_slide(prs, "Additional Analysis Methods", [
        "BCa Bootstrap Confidence Intervals",
        ("Bias-corrected and accelerated", 1),
        ("Handles skewness in IV estimates", 1),
        "Power Analysis",
        ("Required sample size for target power", 1),
        ("Minimum detectable effect (MDE)", 1),
        "Heterogeneity Analysis",
        ("Stratified IV by region (core/margin/infiltrating)", 1),
        ("Stratified by hypoxia (normoxic/mild/severe)", 1),
        ("Cochran's Q test, I² statistic", 1),
    ])

    # ==================== VALIDATION ====================
    add_section_slide(prs, "Validation & Success Criteria")

    add_table_slide(prs, "Success Criteria",
        ["Criterion", "Target", "Measurement"],
        [
            ["IV Accuracy", "< 15% bias", "|β̂_IV - β_true| / β_true"],
            ["CI Coverage", "≥ 90%", "Fraction containing true value"],
            ["Power", "≥ 80%", "At effect size 0.05, n=1000"],
            ["Segregation", "p > 0.05", "KS test vs Binomial(N, 0.5)"],
            ["First-stage F", "> 10", "Staiger & Stock rule"],
            ["E-value", "> 2.0", "Robustness to confounding"],
        ]
    )

    add_content_slide(prs, "Biological Validation Targets", [
        "Avascular spheroid zonation",
        ("Proliferating rim, quiescent zone, necrotic core", 1),
        "Oxygen gradient: hypoxia at ~100-150 μm from vessels",
        "Pseudopalisading necrosis (GBM pathognomonic)",
        "ecDNA heterogeneity dynamics",
        ("Variance = N × p × (1-p) per generation", 1),
        "Angiogenic switch at ~1-2 mm tumor mass",
        "Gompertzian growth curve",
        "Causal recoverability above chance",
    ])

    # ==================== REPOSITORY ====================
    add_section_slide(prs, "Repository Structure")

    add_two_column_slide(prs, "Module Organization",
        "causanta/simulate/",
        [
            "• core.py - Main simulation loop",
            "• cells.py - Cell agents",
            "• behaviors.py - Division, migration, death",
            "• environment.py - Diffusion solver",
            "• ecdna.py - SIV mechanism",
            "• angiogenesis.py - Vessel sprouting",
            "• sweep.py - Parameter sweeps",
        ],
        "causanta/analyze/",
        [
            "• iv.py - 2SLS + diagnostics",
            "• effects.py - Effect estimation",
            "• bootstrap.py - BCa CIs",
            "• power.py - Power analysis",
            "• heterogeneity.py - Stratified IV",
            "• loader.py - Data loading",
            "• discovery.py - PC/GES algorithms",
        ]
    )

    add_content_slide(prs, "Output Files", [
        "cells_t{step}.tsv - Cell states per timestep",
        ("cell_id, parent_id, type, x, y, ecDNA_count, O2_local, ...", 1),
        "environment_t{step}.tsv - Field values",
        ("x_grid, y_grid, O2, glucose, VEGF, lactate, pH, ...", 1),
        "lineage.tsv - Division events",
        ("time_hr, parent_id, parent_ecDNA_before/after, daughter_ecDNA", 1),
        "params.json - Frozen parameter copy",
        "summary.log - Per-step statistics",
    ])

    # ==================== USAGE ====================
    add_section_slide(prs, "Quick Start")

    add_code_slide(prs, "Running Simulations", """
    # Install
    pip install -e .

    # Run default simulation (1mm x 1mm, 168 hours)
    python -m causanta.simulate.core

    # Custom duration
    python -m causanta.simulate.core --hours 500

    # Custom parameters
    python -m causanta.simulate.core my_params.json

    # Comprehensive analysis
    python scripts/run_comprehensive_analysis.py output/run_*/data/

    # Full validation study
    python scripts/run_simulation_study.py --study full --output-dir results/
""")

    add_code_slide(prs, "Generating Figures", """
    from causanta.visualize import (
        apply_nature_style,
        create_main_figure_1,  # Conceptual framework (DAG)
        create_main_figure_2,  # Instrument validation
        create_main_figure_3,  # Causal estimation (IV vs OLS)
        create_main_figure_4,  # Power analysis
        create_main_figure_5,  # Heterogeneity
        create_main_figure_6,  # Robustness (sensitivity)
    )

    apply_nature_style()
    create_main_figure_1(output_dir="figures/")
""")

    # ==================== DOCUMENTATION ====================
    add_content_slide(prs, "Documentation", [
        "docs/theory.md - Complete theoretical framework",
        ("14 sections: domain, cells, environment, behaviors, SIV, ...", 1),
        "docs/manuscript/manuscript.md - Full publication paper",
        "docs/supplementary_materials.md - Statistical methods",
        "docs/tutorial.md - Step-by-step causal inference guide",
        "docs/analysis_guide.md - Methods and interpretation",
        "docs/immune_system.md - Immune interaction modeling",
        "README.md - Collaborator onboarding document",
    ])

    # ==================== SUMMARY ====================
    add_section_slide(prs, "Summary")

    add_content_slide(prs, "Key Takeaways", [
        "CAUSANTA enables causal inference benchmarking with ground truth",
        "ecDNA segregation provides natural randomization (SIV)",
        "6-phase simulation generates realistic tumor tissue",
        "2SLS/IV analysis recovers causal effects despite confounding",
        "Comprehensive diagnostics validate instrument strength",
        "Full documentation for collaborator onboarding",
        "Publication-ready figures with Nature Methods styling",
    ])

    add_title_slide(
        prs,
        "Questions?",
        "docs/theory.md | docs/manuscript/manuscript.md | README.md"
    )

    # Save
    output_path = Path(__file__).parent.parent / "docs" / "CAUSANTA_Overview.pptx"
    prs.save(output_path)
    print(f"Presentation saved to: {output_path}")
    return output_path


if __name__ == "__main__":
    create_presentation()
