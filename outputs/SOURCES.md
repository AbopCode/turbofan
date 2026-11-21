# SOURCES.md

Only sources that were actually opened during this build are listed under "Opened".
Access date for all: 14 November 2026 (code repositories cloned the same day).

## Opened

| # | Source | URL | What it supplied | IEEE entry |
|---|---|---|---|---|
| 1 | NASA AGTF30 repository (commit a1521a76930df92ae5611f3454fc612a8b574cdf, 2020-06-01) | https://github.com/nasa/AGTF30 | Component maps, constants, controller schedules, Simulink wiring (`AGTF30_eng.mdl`), steady-state initial conditions (`setup_InputsSS.m`). Licence: Apache 2.0. | NASA Glenn Research Center, "Advanced Geared Turbofan 30,000 lbf (AGTF30) engine simulation," GitHub repository, 2020. [Online]. Available: https://github.com/nasa/AGTF30 (accessed Nov. 14, 2026). |
| 2 | NASA T-MATS repository (commit ad6e4d5d0d76c229db4eb72ca40ef58d5ddc4014, 2025-01-15) | https://github.com/nasa/T-MATS | C source of every component (compressor, turbine, burner, nozzle, ducts/static calculation, valve, ambient), gas-property tables, interpolation routines, Simulink library blocks (inlet, splitter, gearbox, shaft). Licence: Apache 2.0. | NASA Glenn Research Center, "Toolbox for the Modeling and Analysis of Thermodynamic Systems (T-MATS)," GitHub repository, 2025. [Online]. Available: https://github.com/nasa/T-MATS (accessed Nov. 14, 2026). |
| 3 | J. W. Chapman and J. S. Litt, "Control design for an advanced geared turbofan engine," NASA/TM-2017-219569 | https://ntrs.nasa.gov/citations/20170011210 (PDF: https://ntrs.nasa.gov/api/citations/20170011210/downloads/20170011210.pdf) | Table 1 performance anchors (cruise and max-thrust rows); shaft naming (paper N1/N2/Nf vs code N2/N3/N1); gear ratio; controller description. | J. W. Chapman and J. S. Litt, "Control design for an advanced geared turbofan engine," NASA Glenn Research Center, Cleveland, OH, USA, NASA/TM-2017-219569, Nov. 2017. [Online]. Available: https://ntrs.nasa.gov/citations/20170011210 |
| 4 | M. Pakmehr and T. Yucelen, "Model reference adaptive control of systems with gain scheduled reference models," arXiv:1403.3738 | https://arxiv.org/abs/1403.3738 (PDF: https://arxiv.org/pdf/1403.3738) | The three JetCat SPT5 turboshaft A and B matrices, Sec. 3.3, eqs. (57)-(59), used as the published-example anchor (illustration only). | M. Pakmehr and T. Yucelen, "Model reference adaptive control of systems with gain scheduled reference models," arXiv:1403.3738, Mar. 2014. [Online]. Available: https://arxiv.org/abs/1403.3738 |

## Referred to but NOT opened by this build

These appear in the project brief and in code comments. They were not opened here, so the
team must open them (and fix page numbers/DOIs) before citing.

- J. W. Chapman et al., "Propulsion system simulation using T-MATS," NASA/TM-2014-218410 (the T-MATS user paper). The T-MATS *source code* (#2) was used instead.
- S. M. Jones, W. J. Haller, and M. T. Tong, "An N+3 technology level reference propulsion system," NASA/TM-2017-219501.
- G. H. Golub and C. F. Van Loan, *Matrix Computations*, 4th ed. (Bauer-Fike theorem): standard textbook result; implemented and verified numerically (5,000 random cases in tests/test_mathcore.py), not read from the book here.
- C. J. Daniele and S. M. Krosel, NASA TP-1388 (finite-difference linearization errors).
- N-CMAPSS dataset (NASA Prognostics Data Repository): not used; `./run.sh fetch-data` tries to download it, but its URL could not be verified from the build sandbox.

No Kaggle copy of any dataset is used or cited.
