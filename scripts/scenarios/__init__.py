"""Biological validation scenarios for CAUSANTA.

Each scenario addresses a specific scientific question:

1. margin_invasion.py
   Question: Does ecDNA causally drive invasion at tumor margins?
   Tests IV vs OLS under hypoxia confounding.

2. drug_response.py
   Question: Does ecDNA predict response to EGFR inhibitors?
   Simulates drug effects and tests dose-response relationships.

3. driver_vs_passenger.py
   Question: Can IV distinguish driver from passenger ecDNA?
   Tests sensitivity/specificity of null vs true effect detection.

4. immune_selection.py
   Question: Does immune pressure select for/against high-ecDNA cells?
   Varies immune strength and tracks ecDNA distribution.
"""

from pathlib import Path

SCENARIOS_DIR = Path(__file__).parent

__all__ = ["SCENARIOS_DIR"]
