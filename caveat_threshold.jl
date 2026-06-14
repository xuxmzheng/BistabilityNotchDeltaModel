_t_start = time()   # --- start runtime timer ---
# =====================================================================
#  Caveat figure: operating point pinned at the Hill threshold x = 1.
#  Because H(1)=1/2 for every p, the equilibrium (Nbar,Dbar) is
#  p-independent, so the operating point stays at x = bi*Nbar = 1.
#  The full->QSS error then scales like H'(1)=p/4 and GROWS with p,
#  the opposite of the off-threshold (generic) case.
#
#  Constructed parameter set (dimensionless):
#     ad=2, an=2, bd=7.3, bn=4.6, bi=5,  equilibrium (Nbar,Dbar)=(0.2,0.1)
#  Requires: DifferentialEquations, CairoMakie
#  Run:      julia caveat_threshold.jl  ->  caveat_threshold.png
# =====================================================================
using DifferentialEquations
using CairoMakie

function caveat_threshold()
    # ---- parameters: equilibrium sits exactly at x = bi*Nbar = 1 ----
    ad, an = 2.0, 2.0
    Nb, Db = 0.2, 0.1
    bi = 1.0 / Nb                       # = 5  =>  x = bi*Nb = 1
    bd = (1.5 - ad*Nb*Db)/Nb            # = 7.3
    bn = (0.5 - an*Nb*Db)/Db            # = 4.6
    dI0 = 0.5                           # fixed NICD perturbation: I(0) = bi*Nb + dI0

    Hp(I, p) = 1.0 / (1.0 + I^p)

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

    # e(t) = || (N,D)_full - (N,D)_QSS ||   (QSS starts at the equilibrium and stays)
    function errcurve(p, tau; T = 30.0)
        qsol = solve(ODEProblem(qss!,  [Nb, Db],          (0.0, T), (p = p,)),          Tsit5();  reltol = 1e-10, abstol = 1e-12)
        fsol = solve(ODEProblem(full!, [Nb, Db, bi*Nb+dI0], (0.0, T), (p = p, tau = tau)), Rodas5(); reltol = 1e-9,  abstol = 1e-11)
        ts = range(0, T; length = 2000)
        e  = [sqrt((fsol(t)[1]-qsol(t)[1])^2 + (fsol(t)[2]-qsol(t)[2])^2) for t in ts]
        return collect(ts), e
    end

    fig = CairoMakie.Figure(size = (1300, 560), fontsize = 26)

    # ---------------- panel (a): e(t) for several p (peak grows with p) ----------------
    axa = CairoMakie.Axis(fig[1, 1],
        titlesize = 28, xlabelsize = 26, ylabelsize = 26, xticklabelsize = 22, yticklabelsize = 22,
        title = L"\text{(a)  error vs time at the threshold } x=1",
        xlabel = L"\text{time } t \text{ (dimensionless)}",
        ylabel = L"\Vert(N,D)-(\overline{N},\overline{D})\Vert(t)")
    pp    = [1.0, 2.0, 3.0, 4.0, 5.0]   # p <= 5: threshold equilibrium is stable (p* ~ 5.54)
    stys  = [:solid, :dash, :dot, :dashdot, :dashdotdot]
    mrks  = [:circle, :rect, :diamond, :utriangle, :star5]
    for i in eachindex(pp)
        ts, e = errcurve(pp[i], 0.2; T = 12.0)
        tts = range(0.0, 12.0, length = 16)
        idx = clamp.(round.(Int, collect(tts) ./ 12.0 .* (length(ts) - 1)) .+ 1, 1, length(ts))
        CairoMakie.scatterlines!(axa, collect(tts), e[idx], color = :black, linewidth = 3.5,
            linestyle = stys[i], marker = mrks[i], markersize = 14, label = L"p = %$(pp[i])")
    end
    CairoMakie.xlims!(axa, 0, 6)
    CairoMakie.axislegend(axa, position = (1.0, 0.82), labelsize = 26, patchsize = (58, 22))

    # ---------------- panel (b): sup error vs p, with H'(1)=p/4 prediction ----------------
    ps   = [1.0, 1.5, 2.0, 3.0, 4.0, 5.0]   # p <= 5: equilibrium stable, error is the clean boundary-layer effect
    supe = [maximum(errcurve(p, 0.2; T = 30.0)[2]) for p in ps]
    sens = ps ./ 4.0                                          # H'(1) = p/4
    # single global scale (log-space) so the prediction is not pinned to one point
    pred = sens .* exp(sum(log.(supe) .- log.(sens)) / length(ps))

    axb = CairoMakie.Axis(fig[1, 2],
        titlesize = 28, xlabelsize = 26, ylabelsize = 26, xticklabelsize = 22, yticklabelsize = 22,
        title = L"\text{(b)  sup error grows with } p \text{ at } x=1",
        xlabel = L"\text{Hill coefficient } p", ylabel = L"\mathrm{sup}_t\,\Vert(N,D)-(\overline{N},\overline{D})\Vert(t)")
    CairoMakie.scatterlines!(axb, ps, supe, color = :black, linewidth = 3.5,
                             marker = :circle, markersize = 15, label = L"\text{measured sup error}")
    CairoMakie.lines!(axb, ps, pred, color = :black, linewidth = 3.0,
                      linestyle = :dash, label = L"\text{slope} = 1/4")
    CairoMakie.axislegend(axb, position = :lt, labelsize = 26, patchsize = (58, 22))

    CairoMakie.save("caveat_threshold.png", fig)
    println("[caveat] saved -> ", abspath("caveat_threshold.png"))
    return
end

caveat_threshold()

println("Running time: ", round(time() - _t_start, digits=2), " seconds")
