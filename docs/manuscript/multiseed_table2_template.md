This template is meant to replace the current single-seed point estimates in Table 2 of `docs/manuscript/manuscript.md` (the table in the Robustness Analysis section). Once both batch runs have finished, run `python scripts/aggregate_multiseed_results.py` on your local machine to produce `output/multiseed_full_results.json`, then substitute each `TBD` cell below with the corresponding `mean ± SD` values from that JSON. Update the surrounding prose to reference n=5 reliability bounds rather than single-seed point estimates. When the manuscript has been updated and the draft PR is ready for review, delete this file.

| Scale | Scenario | n_seeds | OLS β (mean ± SD) | IV β (mean ± SD) | OLS bias % (mean ± SD) | F1 (mean ± SD) |
|-------|----------|---------|-------------------|------------------|------------------------|----------------|
| 2mm   | baseline | 5       | TBD               | TBD              | TBD                    | TBD            |
| 2mm   | reduced  | 5       | TBD               | TBD              | TBD                    | TBD            |
| 2mm   | removed  | 5       | TBD               | TBD              | TBD                    | TBD            |
| 6mm   | baseline | 5       | TBD               | TBD              | TBD                    | TBD            |
| 6mm   | reduced  | 5       | TBD               | TBD              | TBD                    | TBD            |
| 6mm   | removed  | 5       | TBD               | TBD              | TBD                    | TBD            |
