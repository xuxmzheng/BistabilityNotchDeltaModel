using DifferentialEquations
using CairoMakie
_t_start = time()   # --- start runtime timer ---

function run_simulation()

    # ---------------- Parameters ----------------
    β_N = 500.0
    β_D = 1000.0
    γ_I = 0.5
    γ = 0.1
    k_T = 0.00005
    k_c = 0.0005
    hill = 1
    Next = 500.0
    Dext = 1110.0
    I_th = 200.0

    N_init = 200.0
    D_init = 200.0
    I_init = 20.0

    H_plus(x)  = (1 + 2*x^hill) / (1 + x^hill)
    H_minus(x) = 1 / (1 + x^hill)

    # -------- Full 3D Model --------
    function f1!(du, u, p, t)
        N, D, I = u
        du[1] = β_N * H_plus(I / I_th) - k_c * N * D - k_T * N * Dext - γ * N
        du[2] = β_D * H_minus(I / I_th) - k_c * N * D - k_T * D * Next - γ * D
        du[3] = k_T * N * Dext - γ_I * I
    end

    # -------- QSS 2D Model --------
    I_expr(N) = k_T * N * Dext / γ_I

    function f2!(du, u, p, t)
        N, D = u
        I = I_expr(N)
        du[1] = β_N * H_plus(I / I_th) - k_c * N * D - k_T * N * Dext - γ * N
        du[2] = β_D * H_minus(I / I_th) - k_c * N * D - k_T * D * Next - γ * D
    end

    tspan = (0.0, 2500.0)

    prob1 = ODEProblem(f1!, [N_init, D_init, I_init], tspan)
    sol1 = solve(prob1, Rodas5(), reltol=1e-8, abstol=1e-8)

    prob2 = ODEProblem(f2!, [N_init, D_init], tspan)
    sol2 = solve(prob2, Rodas5(), reltol=1e-8, abstol=1e-8)

    # nondimensionalization (Sec. 2): Ñ = N/(β_N/γ), D̃ = D/(β_D/γ), Ĩ = I/I_th, t̃ = γ·t
    N0scale = β_N / γ;  D0scale = β_D / γ;  I0scale = I_th
    twin = 2500.0                                  # physical window (hours), only to evaluate the solution
    tdense = range(0, twin, length=400)
    tmark  = range(0, twin, length=15)
    # dimensionless trajectories (functions of physical time t)
    N1(t)=sol1(t)[1]/N0scale; D1(t)=sol1(t)[2]/D0scale; I1(t)=sol1(t)[3]/I0scale
    N2(t)=sol2(t)[1]/N0scale; D2(t)=sol2(t)[2]/D0scale; I2(t)=I_expr(sol2(t)[1])/I0scale
    xd = collect(tdense) .* γ                      # dimensionless time (smooth curves)
    xm = collect(tmark)  .* γ                      # dimensionless time (markers)

    # ---------------- Plot ----------------
    # fig = CairoMakie.Figure(size=(1050,700), fontsize=24, figure_padding=(10,85,10,10))
    # fig = CairoMakie.Figure(size=(1050,700), fontsize=24, figure_padding=(10,10,10,10))
    fig = CairoMakie.Figure(size=(1000,700), fontsize=24, figure_padding=(10,10,10,10))

    ax = CairoMakie.Axis(fig[1,1],
        xlabel = L"\text{time } t \text{ (dimensionless)}",
        ylabel = L"\text{dimensionless concentration}",
        xticks = 0:50:250,
        xticklabelsize = 30,
        yticklabelsize = 30,
        xlabelsize = 42,
        ylabelsize = 42,
        tellwidth = true,
        tellheight = true
    )

    xticks = (0:50:200, string.(0:50:200))
    
    CairoMakie.xlims!(ax, 0, twin*γ*1.05)   # +5% so the last tick (250) is not clipped at the edge
    CairoMakie.ylims!(ax, 0, nothing)
    
    blk = :black
    # colorblind-safe monochrome: all curves black,
    #   Full = solid line + filled marker,  QSS = dashed line + open marker,
    #   variable = marker shape (circle = N, square = D, triangle = I)

    # ---- Full model (solid line, filled markers) ----
    CairoMakie.lines!(ax, xd, N1.(tdense), color=blk, linewidth=4.5, linestyle=:solid)
    CairoMakie.lines!(ax, xd, D1.(tdense), color=blk, linewidth=4.5, linestyle=:solid)
    CairoMakie.lines!(ax, xd, I1.(tdense), color=blk, linewidth=4.5, linestyle=:solid)
    CairoMakie.scatter!(ax, xm, N1.(tmark), color=blk, marker=:circle,    markersize=28)
    CairoMakie.scatter!(ax, xm, D1.(tmark), color=blk, marker=:rect,      markersize=27)
    CairoMakie.scatter!(ax, xm, I1.(tmark), color=blk, marker=:utriangle, markersize=29)

    # ---- QSS model (dashed line, open markers) ----
    CairoMakie.lines!(ax, xd, N2.(tdense), color=blk, linewidth=4.5, linestyle=:dash)
    CairoMakie.lines!(ax, xd, D2.(tdense), color=blk, linewidth=4.5, linestyle=:dash)
    CairoMakie.lines!(ax, xd, I2.(tdense), color=blk, linewidth=4.5, linestyle=:dash)
    CairoMakie.scatter!(ax, xm, N2.(tmark), color=:white, strokecolor=blk, strokewidth=1.5, marker=:circle,    markersize=28)
    CairoMakie.scatter!(ax, xm, D2.(tmark), color=:white, strokecolor=blk, strokewidth=1.5, marker=:rect,      markersize=27)
    CairoMakie.scatter!(ax, xm, I2.(tmark), color=:white, strokecolor=blk, strokewidth=1.5, marker=:utriangle, markersize=29)

    # ---- legend: 2 rows (orientation=:horizontal, nbanks=2) -> row 1 = Full (N,D,I), row 2 = QSS (N,D,I) ----
    LE(ls) = CairoMakie.LineElement(color=blk, linewidth=4.5, linestyle=ls)
    MEf(m) = CairoMakie.MarkerElement(color=blk, marker=m, markersize=27)
    MEo(m) = CairoMakie.MarkerElement(color=:white, strokecolor=blk, strokewidth=1.5, marker=m, markersize=27)
    # interleaved so a horizontal 2-bank legend gives row1 = Full(N,D,I), row2 = QSS(N,D,I)
    elems = [[LE(:solid), MEf(:circle)],    [LE(:dash), MEo(:circle)],
             [LE(:solid), MEf(:rect)],      [LE(:dash), MEo(:rect)],
             [LE(:solid), MEf(:utriangle)], [LE(:dash), MEo(:utriangle)]]
    labels = [L"N \text{ (Full)}", L"N \text{ (QSS)}", L"D \text{ (Full)}", L"D \text{ (QSS)}", L"I \text{ (Full)}", L"I \text{ (QSS)}"]
    CairoMakie.axislegend(ax, elems, labels, position = (1.0, 0.72), labelsize = 40, patchsize = (62, 28), nbanks = 2, orientation = :horizontal)

    CairoMakie.save("compare_full_QSS_p=1.png", fig)
    println("[colorblind v] saved black/marker figure -> ", abspath("compare_full_QSS_p=1.png"))
    display(fig)

end

run_simulation()

println("Running time: ", round(time() - _t_start, digits=2), " seconds")
