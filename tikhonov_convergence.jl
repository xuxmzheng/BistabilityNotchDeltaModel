# =====================================================================
#  Tikhonov / singular-perturbation convergence test
#  Full Notch-Delta-NICD system (non1)  -->  QSS system (nonC1)  as tau -> 0
#
#  Confirms numerically:
#   (i)  LONG-TIME convergence: full and QSS solutions reach the SAME
#        equilibrium from the same slow initial data (both basins).
#   (ii) sup_t || (N,D)_full(t) - (N,D)_QSS(t) ||  =  O(tau)   (log-log slope ~ 1).
#
#  Model (dimensionless, slow time t; see paper Eqs. (non1)/(nonC1)):
#     dN/dt      = 2 - 1/(1+I^p) - ad*N*D - bd*N
#     dD/dt      =     1/(1+I^p) - an*N*D - bn*D
#     tau*dI/dt  = bi*N - I                 (tau = gamma/gamma_I)
#  QSS limit (tau->0): I = bi*N substituted into the (N,D) equations.
#
#  Fixed dimensionless parameters (from paper Table 1):
#     ad = kC*betaD/gamma^2 = 50,  an = kC*betaN/gamma^2 = 25,  p = 2
#     bd = kT*Dext/gamma + 1,  bn = kT*Next/gamma + 1,  bi = betaN*kT*Dext/(I0*gammaI*gamma)
#  Regime: Bistable  Dext=800, Next=200  (3 equilibria: two stable, one saddle).
#
#  Requires:  DifferentialEquations, Plots, LinearAlgebra, Printf, Statistics
#  Run:       julia tikhonov_convergence.jl     ->  tikhonov_convergence_julia.png
# =====================================================================

using DifferentialEquations
using LinearAlgebra
using Printf
using Statistics
using Plots
using Plots.PlotMeasures      # provides mm for margins (ships with Plots; no extra pkg)
using LaTeXStrings            # L"..." math labels  (add with: ] add LaTeXStrings)

# matplotlib backend: its legend line sample is long enough to show dash patterns
# (GR draws too short a sample, so :dash and :dashdot both look solid in the legend).
import PythonPlot              # add with: ] add PythonPlot   (needs matplotlib)
_t_start = time()   # --- start runtime timer ---
pythonplot()
PythonPlot.rc("legend", handlelength = 4.0)   # lengthen the legend line sample

# larger legends and axis (guide) labels for all panels;
# per-subplot margins (set here so EACH panel inherits them) keep the titles /
# axis labels of adjacent panels from colliding in the 2x3 grid
Plots.default(guidefontsize=22, legendfontsize=18, tickfontsize=12, titlefontsize=19,
        framestyle=:box, grid=true, gridalpha=0.18, gridstyle=:dot,
        foreground_color_legend=nothing,
        left_margin=6mm, right_margin=4mm, top_margin=5mm, bottom_margin=9mm)

# colorblind-safe (Okabe-Ito) colors + distinct line styles for the three tau curves
const FULL_STY = [:dash, :dot, :dashdot]
const FULL_COL = [:black, :black, :black]   # black-and-white: curves distinguished by line style only

const AD = 50.0
const AN = 25.0
const PP = 2.0

# ---- parameter map from external Delta/Notch -------------------------
function make_params(Dext, Next)
    bd = 5e-4*Dext + 1.0
    bn = 5e-4*Next + 1.0
    bi = 0.0025*Dext
    return (ad=AD, an=AN, p=PP, bd=bd, bn=bn, bi=bi)
end

Hp(x, p) = 1.0 / (1.0 + x^p)

# ---- vector fields ---------------------------------------------------
function full!(du, u, P, t)        # P = (params..., tau)
    N, D, I = u
    du[1] = 2 - Hp(I, P.p) - P.ad*N*D - P.bd*N
    du[2] =     Hp(I, P.p) - P.an*N*D - P.bn*D
    du[3] = (P.bi*N - I) / P.tau
    return nothing
end

function qss!(du, u, P, t)
    N, D = u
    h = Hp(P.bi*N, P.p)
    du[1] = 2 - h - P.ad*N*D - P.bd*N
    du[2] =     h - P.an*N*D - P.bn*D
    return nothing
end

# ---- one convergence study for a given regime & initial condition ----
function study(Dext, Next, label, IC2, I0;
               T = 25.0,
               taus_plot = [0.2, 0.1, 0.025],
               taus_err  = [0.4, 0.2, 0.1, 0.05, 0.025, 0.0125, 0.00625, 0.003125])
    base = make_params(Dext, Next)

    # reference QSS solution (dense)
    qprob = ODEProblem(qss!, [IC2[1], IC2[2]], (0.0, T), base)
    qsol  = solve(qprob, Tsit5(); reltol=1e-10, abstol=1e-12)
    Nbar(t) = qsol(t)[1];  Dbar(t) = qsol(t)[2]

    # --- sup-norm error vs tau ---
    errs = Float64[]
    for tau in taus_err
        P = merge(base, (tau=tau,))
        fprob = ODEProblem(full!, [IC2[1], IC2[2], I0], (0.0, T), P)
        fsol  = solve(fprob, Rodas5P(); reltol=1e-9, abstol=1e-11)   # stiff solver (tau small)
        ts = range(0, T; length=4000)
        e  = maximum(@. sqrt((getindex.(fsol.(ts),1) - Nbar(ts))^2 +
                             (getindex.(fsol.(ts),2) - Dbar(ts))^2))
        push!(errs, e)
    end
    A     = hcat(log.(taus_err[end-3:end]), ones(4))
    slope = (A \ log.(errs[end-3:end]))[1]
    @printf("%-24s : log-log slope (small-tau) = %.3f\n", label, slope)

    # --- column 1: N(t) full at several tau vs QSS ---
    pN = Plots.plot(title=label, xlabel=L"time\ t\ (dimensionless)", ylabel=L"N(t),\ dimensionless", legend=:best)
    tt = range(0, T; length=600)
    Plots.plot!(pN, tt, Nbar.(tt), lw=3, ls=:solid, color=:black, label=L"QSS\ \ \overline{N}(t)")
    for (i, tau) in enumerate(taus_plot)
        P = merge(base, (tau=tau,))
        fprob = ODEProblem(full!, [IC2[1], IC2[2], I0], (0.0, T), P)
        fsol  = solve(fprob, Rodas5P(); reltol=1e-9, abstol=1e-11)
        Plots.plot!(pN, tt, [fsol(t)[1] for t in tt], lw=2.4, ls=FULL_STY[i], color=FULL_COL[i],
              label=L"full\ \tau=%$tau")
    end

    # --- column 2: boundary layer in I(t) ---
    tw = 2.2
    pI = Plots.plot(title=L"boundary\ layer\ in\ I(t)", xlabel=L"time\ t\ (dimensionless)", ylabel=L"I(t),\ dimensionless", legend=:best)
    tw_grid = range(0, tw; length=600)
    Plots.plot!(pI, tw_grid, base.bi .* Nbar.(tw_grid), lw=2.8, color=:black,
          label=L"critical\ manifold\ \beta_i\, \overline{N}(t)")
    for (i, tau) in enumerate(taus_plot)
        P = merge(base, (tau=tau,))
        fprob = ODEProblem(full!, [IC2[1], IC2[2], I0], (0.0, tw), P)
        fsol  = solve(fprob, Rodas5P(); reltol=1e-9, abstol=1e-11)
        Plots.plot!(pI, tw_grid, [fsol(t)[3] for t in tw_grid], lw=2.4, ls=FULL_STY[i], color=FULL_COL[i],
              label=L"full\ I(t),\ \tau=%$tau")
    end

    # --- column 3: error scaling (log-log) ---
    pE = Plots.plot(xscale=:log10, yscale=:log10, legend=:topleft, legendfontsize=16,
              xticks=([0.01, 0.1, 0.2], ["0.01", "0.1", "0.2"]),
              xlabel=L"\tau = \gamma/\gamma_I", ylabel=L"\mathrm{sup}_t\,\Vert(N,D)-(\overline{N},\overline{D})\Vert(t)",
              title=L"error\ scaling\ (slope\approx %$(round(slope,digits=2)))", yguidefontsize=13)
    # dummy NaN point -> controls ONLY the (small) "sup error" symbol in the legend
    Plots.scatter!(pE, [NaN], [NaN], markersize=1.5, markerstrokewidth=0.2, color=:black, label="sup error")
    # actual data points (markers only, no connecting line), no legend entry
    Plots.scatter!(pE, taus_err, errs, markersize=7.0, markerstrokewidth=0.5, color=:black, label="")
    ref = errs[end] / taus_err[end] .* taus_err          # slope-1 guide
    Plots.plot!(pE, taus_err, ref, ls=:dot, color=:black, lw=3.0, label=L"slope\ 1\ \ (O(\tau))")
    Plots.vline!(pE, [0.2], ls=:dash, color=:gray, lw=2.5, label=L"\tau=0.2")

    return pN, pI, pE
end

# ============================ run ====================================
pN1, pI1, pE1 = study(800, 200, L"Bistable:\ low{-}N\ basin",  (0.05, 0.5), 3.0)
pN2, pI2, pE2 = study(800, 200, L"Bistable:\ high{-}N\ basin", (1.2, 0.05), 0.0)

# high-N basin: move the N(t) legend to the bottom right
Plots.plot!(pN2, legend=:bottomright, legendfontsize=14)

# 2nd column (I(t)): smaller legends -- top row top-right, bottom row bottom-right
Plots.plot!(pI1, legend=:topright,    legendfontsize=15)
Plots.plot!(pI2, legend=:bottomright, legendfontsize=15)

# layout: 2 rows x 3 cols -- row 1 = low-N basin, row 2 = high-N basin
fig = Plots.plot(pN1, pI1, pE1, pN2, pI2, pE2; layout=(2,3), size=(1600, 900), dpi=200)
Plots.savefig(fig, "tikhonov_convergence.png")
println("\nSaved figure: tikhonov_convergence.png")

println("Running time: ", round(time() - _t_start, digits=2), " seconds")
