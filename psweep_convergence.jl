# =====================================================================
#  Fig. 10 : Impact of the Hill coefficient p on full -> QSS convergence
#  (dimensionless QSS test, Fig.8 parameter set: Dext=1110, Next=500)
#
#  Panel (a): full-QSS slow error vs time at tau=0.2 for several p.
#             COLORBLIND-SAFE: all curves BLACK, distinguished by line style.
#  Panel (b): peak error (left axis) and coincidence time (right axis) vs p.
#  Panel (c): O(tau) prefactor C(p) vs the Hill sensitivity prediction.
#
#  Requires: DifferentialEquations, CairoMakie
#  Run:      julia psweep_convergence.jl   ->  psweep_convergence.png/.pdf
# =====================================================================
using DifferentialEquations
using CairoMakie
_t_start = time()   # --- start runtime timer ---

function psweep_convergence()
    # ---- dimensionless parameters (Dext=1110, Next=500) ----
    ad = 50.0; an = 25.0
    bd = 5e-4*1110 + 1.0      # 1.555
    bn = 5e-4*500  + 1.0      # 1.25
    bi = 0.0025*1110          # 2.775
    N0, D0, I0 = 0.04, 0.02, 0.1

    Hp(x, p) = 1.0 / (1.0 + x^p)

    function full!(du, u, par, t)
        N, D, I = u
        du[1] = 2 - Hp(I, par.p) - ad*N*D - bd*N
        du[2] =     Hp(I, par.p) - an*N*D - bn*D
        du[3] = (bi*N - I) / par.tau
    end
    function qss!(du, u, par, t)
        N, D = u
        h = Hp(bi*N, par.p)
        du[1] = 2 - h - ad*N*D - bd*N
        du[2] =     h - an*N*D - bn*D
    end

    # e(t) = || (N,D)_full - (N,D)_QSS ||  for given p, tau
    function errcurve(p, tau; T = 200.0)
        qsol = solve(ODEProblem(qss!,  [N0, D0],     (0.0, T), (p = p,)),          Tsit5();  reltol = 1e-10, abstol = 1e-12)
        fsol = solve(ODEProblem(full!, [N0, D0, I0], (0.0, T), (p = p, tau = tau)), Rodas5(); reltol = 1e-9,  abstol = 1e-11)
        ts = range(0, T; length = 2000)
        e  = [sqrt((fsol(t)[1]-qsol(t)[1])^2 + (fsol(t)[2]-qsol(t)[2])^2) for t in ts]
        return collect(ts), e
    end

    fig = CairoMakie.Figure(size = (1850, 640), fontsize = 28)

    # ---------------- panel (a): error(t), several p, tau=0.2 ----------------
    #   colorblind: black curves, distinct line styles
    axa = CairoMakie.Axis(fig[1, 1], yscale = log10,
        titlesize = 30, xlabelsize = 28, ylabelsize = 28, xticklabelsize = 24, yticklabelsize = 24,
        title = L"\text{(a)  full-QSS error vs time } (\tau = 0.2)",
        xlabel = L"\text{time } t \text{ (dimensionless)}", ylabel = L"\Vert(N,D)-(\overline{N},\overline{D})\Vert(t)")
    pp      = [1.0, 2.0, 3.0, 4.0, 6.0]
    # pp      = [0.1, 0.2, 0.5, 0.7, 0.9]
    stys    = [:solid, :dash, :dot, :dashdot, :dashdotdot]   # 5 named styles; 6th reuses :solid (marker disambiguates)
    markers = [:circle, :rect, :diamond, :utriangle, :dtriangle]
    # Draw each curve as ONE labeled scatterlines! (line + markers); the legend is then
    # auto-generated from the plotted objects, so it always matches the curves.
    for i in eachindex(pp)
        ts, e = errcurve(pp[i], 0.2)                       # ts = range(0, 40, length = 2000)
        tts = range(0.0, 5.0, length = 16)                # coarse grid: smooth-enough line + visible markers
        idx = clamp.(round.(Int, collect(tts) ./ 40.0 .* (length(ts) - 1)) .+ 1, 1, length(ts))
        ee  = e[idx] .+ 1e-12
        CairoMakie.scatterlines!(axa, collect(tts), ee, color = :black, linewidth = 4.0,
            linestyle = stys[i], marker = markers[i], markersize = 18, label = L"p = %$(pp[i])")
    end
    CairoMakie.xlims!(axa, 0, 4); CairoMakie.ylims!(axa, 1e-10, 1)
    CairoMakie.axislegend(axa, position = (1.0, 0.82), labelsize = 28, patchsize = (100, 20), patchlabelgap = 8)

    # ---------------- panel (b): peak error & coincidence time vs p ----------------
    ps = [0.1, 0.5, 0.7, 1.0, 1.3, 1.5, 1.8, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0]
    peak = Float64[]; tco = Float64[]
    for p in ps
        ts, e = errcurve(p, 0.2)
        push!(peak, maximum(e))
        idx = findlast(>(1e-4), e)
        push!(tco, idx === nothing ? 0.0 : ts[idx])
    end
    axb = CairoMakie.Axis(fig[1, 2], yscale = log10,
        titlesize = 30, xlabelsize = 28, ylabelsize = 28, xticklabelsize = 24, yticklabelsize = 24,
        title = L"\text{(b)  sup error and coincidence time}",
        xlabel = L"\text{Hill coefficient } p", ylabel = L"\mathrm{sup}_t\,\Vert(N,D)-(\overline{N},\overline{D})\Vert(t)",
        xticks = [1, 2, 4, 6, 8])
    CairoMakie.scatterlines!(axb, ps, peak, color = :black, linewidth = 4.0,
                             marker = :circle, markersize = 20, label = "sup error")
    # p = 1 marks the maximum of peak error vs p: increasing for p<1, decreasing for p>1
    # (the dotted line meets the x-axis exactly at the tick p = 1)
    CairoMakie.vlines!(axb, [1.0], color = :gray, linestyle = :dot, linewidth = 3)
    axb2 = CairoMakie.Axis(fig[1, 2], yaxisposition = :right, ylabel = L"\text{coincidence time (dimensionless)}",
        ylabelsize = 28, yticklabelsize = 24)
    CairoMakie.hidespines!(axb2); CairoMakie.hidexdecorations!(axb2)
    CairoMakie.linkxaxes!(axb, axb2)
    CairoMakie.scatterlines!(axb2, ps, tco, color = :black, linewidth = 3.5,
                             linestyle = :dash, marker = :rect, markersize = 19, label = "coincidence time")
    # legend with BOTH series (peak error + coincidence time), top-right, larger
    LEb(ls) = CairoMakie.LineElement(color = :black, linewidth = 4.0, linestyle = ls)
    MEb(m)  = CairoMakie.MarkerElement(color = :black, marker = m, markersize = 20)
    CairoMakie.axislegend(axb,
        [[LEb(:solid), MEb(:circle)], [LEb(:dash), MEb(:rect)]],
        [L"\text{sup error}", L"\text{coincidence time}"],
        position = :rt, labelsize = 28, patchsize = (60, 18))

    # ---------------- panel (c): O(tau) prefactor C(p) ----------------
    ps2 = [0.1, 0.5, 0.7, 1.0, 1.3, 1.5, 1.8, 2.0, 3.0, 4.0, 5.0, 6.0, 8.0]
    taus = [0.2, 0.1, 0.05, 0.025, 0.0125]
    Cps  = Float64[]
    for p in ps2
        # Use only the TWO SMALLEST tau (deepest in the asymptotic O(tau) regime).
        # e(tau) = C*tau + O(tau^2)  =>  ratio r = e/tau = C + B*tau.
        # The line through the two ratio points hits tau=0 at the prefactor C
        # (2-point Richardson extrapolation; eliminates the leading O(tau) term).
        t1, t2 = taus[end-1], taus[end]                   # 0.025, 0.0125
        e1 = maximum(errcurve(p, t1)[2])                  # peak slow error at t1
        e2 = maximum(errcurve(p, t2)[2])                  # peak slow error at t2
        r1, r2 = e1 / t1, e2 / t2                         # ratios -> C as tau->0
        B = (r1 - r2) / (t1 - t2)                         # slope of r vs tau
        C = r2 - B * t2                                   # intercept = tau->0 prefactor
        push!(Cps, C)
    end
    xe   = bi * 0.046                                   # operating point beta_i * N_eq ≈ 0.12
    sens = [p * xe^(p-1) / (1 + xe^p)^2 for p in ps2]
    # Scale the Hill-sensitivity prediction to the measured C(p) by ONE global
    # factor (log-space least squares = geometric-mean alignment), so the dashed
    # curve depends only on the SHAPE sens(p), not on which p is listed first.
    logscale = sum(log.(Cps) .- log.(sens)) / length(ps2)
    pred = sens .* exp(logscale)
    axc = CairoMakie.Axis(fig[1, 3], yscale = log10,
        titlesize = 30, xlabelsize = 28, ylabelsize = 28, xticklabelsize = 24, yticklabelsize = 24,
        title = L"(c)  O(\tau) \text{ prefactor } C(p)",
        xlabel = L"\text{Hill coefficient } p", ylabel = L"C(p)")
    CairoMakie.scatterlines!(axc, ps2, Cps, color = :black, linewidth = 4.0,
                             marker = :circle, markersize = 20, label = L"\text{measured } C(p) = \text{error}/\tau")
    CairoMakie.lines!(axc, ps2, pred, color = :black, linewidth = 3.5,
                      linestyle = :dash, label = L"\text{Hill sensitivity} (\beta_i \overline{N} \approx 0.12)")
    CairoMakie.axislegend(axc, position = :lb, labelsize = 28)

    CairoMakie.save("psweep_convergence.png", fig)
    println("[psweep] saved -> ", abspath("psweep_convergence.png"))
    return
end

psweep_convergence()
println("Running time: ", round(time() - _t_start, digits=2), " seconds")
