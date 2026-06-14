# Bistability Analysis of a Hybrid Notch–Delta Signaling Model

Julia codes for the paper *"Bistability Analysis of a Hybrid Notch-Delta Signaling Model in a Single Cell Interacting with Fixed External Proteins"*.

**Authors:** Yuxin Liang, Ling Xue, Kun Zhao, Xiaoming Zheng

Each script is self-contained and writes its figure to a PNG of the same base name. Run a script from the Julia REPL with `include("script.jl")`. The first run in a session is slow because of Julia's just-in-time compilation; subsequent runs are fast. Every script prints its running time when it finishes.

**Requirements:** `DifferentialEquations`, `CairoMakie`, `LaTeXStrings`, `Roots`. The script `tikhonov_convergence.jl` additionally uses `Plots` and `PythonPlot` (matplotlib backend).

## Figures and Generating Scripts

| Fig./Place | Image file | Generating script |
|------------|------------|-------------------|
| Remark 1 (five-root case) | — | `try_5_roots_sweep.jl` |
| Fig 1 | `Notch_Delta_signaling.png` | — (schematic) |
| Fig 2(a) | `bifurcation_p1_Next500_smallwindow.png` | `aaa_p1_Next500_smallwindow.jl` |
| Fig 2(b) | `bifurcation_p2_Next500_smallwindow.png` | `aaa_p2_Next500_smallwindow.jl` |
| Fig 2(c) | `bifurcation_p2_Next12000_smallwindow.png` | `aaa_p2_Next12000_smallwindow.jl` |
| Fig 2(d) | `bifurcation_loglog_p=1_Next=500.png` | `aaa_p1_Next500.jl` |
| Fig 2(e) | `bifurcation_loglog_p=2_Next=500.png` | `aaa_p2_Next500.jl` |
| Fig 2(f) | `bifurcation_loglog_p=2_Next=12000.png` | `aaa_p2_Next12000.jl` |
| Fig 3(a) | `phase_plane_p=1.png` | `phase_plane_p=1.jl` |
| Fig 3(b) | `phase_plane_p=2.png` | `phase_plane_p=2.jl` |
| Fig 4(a) | `bifurcation_p7.7_Next500_smallwindow.png` | `aaa_p7.7_Next500_smallwindow.jl` |
| Fig 4(b) | `bifurcation_p20_Next500_smallwindow.png` | `aaa_p20_Next500_smallwindow.jl` |
| Fig 5(a) | `pheno_p1_BW.png` | `pheno_p1_BW.jl` |
| Fig 5(b) | `pheno_p2_BW.png` | `pheno_p2_BW.jl` |
| Fig 5(c) | `pheno_p4_BW.png` | `pheno_p4_BW.jl` |
| Fig 5(d) | `pheno_p2_special_BW.png` | `pheno_p2_special_BW.jl` |
| Fig 5(e) | `bifurcation_special_p2_lower_BW.png` | `aaa_p2_special_lower_smallwindow_BW.jl` |
| Fig 5(f) | `bifurcation_special_p2_higher_BW.png` | `aaa_p2_special_higher_smallwindow_BW.jl` |
| Fig 6 | `phenotype_diagram.png` | — (schematic) |
| Fig 7(a) | `phenotypebifur_SS05_BW.png` | `phenotypebifur_SS05_BW.jl` |
| Fig 7(b) | `bifur_SS05.png` | `bifur_SS05.jl` |
| Fig 8 | `tikhonov_convergence.png` | `tikhonov_convergence.jl` |
| Fig 9(a) | `compare_full_QSS_p=1.png` | `compare_full_QSS_p=1.jl` |
| Fig 9(b) | `compare_full_QSS_p=2.png` | `compare_full_QSS_p=2.jl` |
| Fig 10 | `psweep_convergence.png` | `psweep_convergence.jl` |
| Fig 11 | `caveat_threshold.png` | `caveat_threshold.jl` |
| Parameter-sweep tables | — | `phenotype_sweep_strategy1.jl`, `phenotype_sweep_strategy2.jl` |
