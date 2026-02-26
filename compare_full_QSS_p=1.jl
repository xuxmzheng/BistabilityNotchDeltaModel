using DifferentialEquations
using CairoMakie

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
    I_init = 0.0

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

    t1 = sol1.t
    t2 = sol2.t
    I2_vals = I_expr.(sol2[1, :])

    # ---------------- Plot ----------------
    fig = CairoMakie.Figure(size = (1000, 700), fontsize = 24)

    ax = CairoMakie.Axis(fig[1,1],
        xlabel = "Time (hours)",
        ylabel = "Molecules",
        xticklabelsize = 30,
        yticklabelsize = 30,
        xlabelsize = 42,
        ylabelsize = 42
    )

    CairoMakie.xlims!(ax, 0, 2550)
    CairoMakie.ylims!(ax, 0, 4500)

    # ---- Full model (solid) ----
    CairoMakie.lines!(ax, t1, sol1[1,:], color=:red,   linewidth=3, label="N (Full)")
    CairoMakie.lines!(ax, t1, sol1[2,:], color=:green, linewidth=3, label="D (Full)")
    CairoMakie.lines!(ax, t1, sol1[3,:], color=:blue,  linewidth=3, label="I (Full)")

    # ---- QSS model (dashed) ----
    CairoMakie.lines!(ax, t2, sol2[1,:], color=:red,   linewidth=3, linestyle=:dash, label="N (QSS)")
    CairoMakie.lines!(ax, t2, sol2[2,:], color=:green, linewidth=3, linestyle=:dash, label="D (QSS)")
    CairoMakie.lines!(ax, t2, I2_vals,  color=:blue,  linewidth=3, linestyle=:dash, label="I (QSS)")

    CairoMakie.axislegend(
        ax,
        position = :rt,
        labelsize = 40,
        patchsize = (40, 20)
    )
    CairoMakie.save("compare_full_QSS_p=1.png", fig)
    display(fig)

end

run_simulation()