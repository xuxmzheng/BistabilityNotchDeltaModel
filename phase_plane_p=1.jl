using DifferentialEquations
using LinearAlgebra
import Plots

Plots.gr()   # Force GR backend

function run_phase_diagram_with_arrows()

    # -------------------------------
    # Parameters
    # -------------------------------
    β_N = 500.0
    β_D = 1000.0
    γ_I = 0.5
    γ = 0.1
    k_T = 0.00005
    k_c = 0.0005
    p = 1
    Dext = 1000.0
    Next = 500.0
    I_threshold = 200.0

    H_plus(x)  = (1 + 2*x^p) / (1 + x^p)
    H_minus(x) = 1 / (1 + x^p)

    # -------------------------------
    # ODE system
    # -------------------------------
    function odefun!(du, u, p, t)
        N, D, I = u
        du[1] = β_N * H_plus(I / I_threshold) - k_c * N * D - k_T * N * Dext - γ * N
        du[2] = β_D * H_minus(I / I_threshold) - k_c * N * D - k_T * D * Next - γ * D
        du[3] = k_T * N * Dext - γ_I * I
    end

    # -------------------------------
    # Initial conditions
    # -------------------------------
    init_conds = [
        [500,   0,    200],
        [1500,  5000, 100],
        [0,     2500, 110],
        [0,     1000, 200],
        [6000,  5000, 200],
        [6000,  3000, 0],
        [1000,  0,    320],
        [3000,  5000, 109],
        [2000,  0,    250],
        [5000,  5000, 80],
        [6000,  2000, 500],
        [3500,  0,    300]
    ]

    colors = [:blue, :blue, :blue, :blue,
              :red,  :red,  :red,
              :blue, :red,
              :blue, :red,  :red]

    tspan = (1.0, 500000.0)

    # -------------------------------
    # Plot setup
    # -------------------------------
    plt = Plots.plot(
        xlabel = "N",
        ylabel = "D",
        xlims = (0, 6300),
        ylims = (0, 5000),
        legend = (0.55, 0.85),
        legendfontsize = 16,
        legend_background_color = :white,
        tickfontsize = 16
    )

    # -------------------------------
    # Solve & plot trajectories
    # -------------------------------
    for (i, u0) in enumerate(init_conds)

        prob = ODEProblem(odefun!, u0, tspan)
        sol  = solve(prob, Rodas5(), reltol=1e-8, abstol=1e-8)

        N_vals = getindex.(sol.u, 1)
        D_vals = getindex.(sol.u, 2)

        label_text = ""
        if i == 1
            label_text = "Receiver basin"
        elseif i == 5
            label_text = "Sender basin"
        end

        Plots.plot!(plt, N_vals, D_vals,
            color = colors[i],
            linewidth = 1.0,
            label = label_text
        )

        # ------------------------------------------------
        # Uniform arc-length based arrow placement
        # ------------------------------------------------
        arc = zeros(length(N_vals))
        for k in 2:length(N_vals)
            arc[k] = arc[k-1] + hypot(
                N_vals[k] - N_vals[k-1],
                D_vals[k] - D_vals[k-1]
            )
        end

        total_length = arc[end]
        n_arrows = 3
        target_positions = range(0, total_length,
                                 length = n_arrows + 2)[2:end-1]

        du = zeros(3)

        for s in target_positions

            idx = findfirst(x -> x ≥ s, arc)
            idx === nothing && continue

            u_arrow = sol.u[idx]
            odefun!(du, u_arrow, nothing, sol.t[idx])

            dN, dD = du[1], du[2]

            # Constant arrow length
            scale = 150.0 / (hypot(dN, dD) + 1e-12)

            Plots.quiver!(plt,
                [u_arrow[1]],
                [u_arrow[2]],
                quiver = ([dN * scale], [dD * scale]),
                color = colors[i],
                linewidth = 1.0,
                label = ""
            )
        end
    end

    # -------------------------------
    # Equilibria
    # -------------------------------
    E1 = [447.623779418578, 2342.57660264528]
    E2 = [1530.0, 600.0]
    E3 = [2448.16352009529, 333.281327303735]

    Plots.scatter!(plt, [E1[1]], [E1[2]],
        markersize=8, color=:blue, label="E₁")

    Plots.scatter!(plt, [E2[1]], [E2[2]],
        markersize=8, color=:black, label="E₂")

    Plots.scatter!(plt, [E3[1]], [E3[2]],
        markersize=8, color=:red, label="E₃")

    Plots.annotate!(plt, 630, 2700, Plots.text("E₁", 20, :blue))
    Plots.annotate!(plt, 1600, 900, Plots.text("E₂", 20, :black))
    Plots.annotate!(plt, 2500, 700, Plots.text("E₃", 20, :red))

    Plots.display(plt)
    Plots.savefig(plt, "phase_plane_p=1.png")

    return plt
end

run_phase_diagram_with_arrows()