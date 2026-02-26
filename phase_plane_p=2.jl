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
    p = 2
    Dext = 1000.0
    Next = 500.0
    I_threshold = 200.0

    H_plus(x)  = (1 + 2*x^p) / (1 + x^p)
    H_minus(x) = 1 / (1 + x^p)

    # -------------------------------
    # ODE System
    # -------------------------------
    function odefun!(du, u, p, t)
        N, D, I = u
        du[1] = β_N * H_plus(I / I_threshold) - k_c * N * D - k_T * N * Dext - γ * N
        du[2] = β_D * H_minus(I / I_threshold) - k_c * N * D - k_T * D * Next - γ * D
        du[3] = k_T * N * Dext - γ_I * I
    end

    # -------------------------------
    # Initial Conditions
    # -------------------------------
    init_conds = [
        [500,   0,    400],
        [1000,  0,    500],
        [4000,  5000, 0],
        [0,     200,  200],
        [500,   5000, 100],
        [6000,  5000, 100],
        [0,     5000, 160],
        [6000,  4000, 100],
        [6000,  500,  600],
        [6000,  0,    400],
        [5000,  0,    400],
        [3000,  0,    400]
    ]

    colors = [:blue, :red, :blue, :blue, :blue, :blue,
              :blue, :red, :red, :red, :red, :red]

    tspan_list = [(1.0, i==1 ? 50000.0 : 500000.0)
                  for i in 1:length(init_conds)]

    # -------------------------------
    # Plot Setup
    # -------------------------------
    plt = Plots.plot(
        xlabel = "N",
        ylabel = "D",
        xlims = (0, 6300),
        ylims = (0, 5000),
        legend = (0.55, 0.85),
        legendfontsize = 16,
        tickfontsize = 16
    )

    # -------------------------------
    # Solve and Plot
    # -------------------------------
    for (i, u0) in enumerate(init_conds)

        prob = ODEProblem(odefun!, u0, tspan_list[i])
        sol  = solve(prob, Rodas5(), reltol=1e-8, abstol=1e-8)

        N_vals = getindex.(sol.u, 1)
        D_vals = getindex.(sol.u, 2)

        label_text = ""
        if i == 1
            label_text = "Receiver basin"
        elseif i == 2
            label_text = "Sender basin"
        end

        Plots.plot!(plt, N_vals, D_vals,
                    color = colors[i],
                    linewidth = 1.0,
                    label = label_text)

        # ------------------------------------------------
        # Uniform Arc-Length Arrow Placement
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

            # Constant arrow size
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
    E1 = [229.204151417523, 4119.48408364746]
    E2 = [2000.0, 400.0]
    E3 = [5550.54030894494, 39.6219875911302]

    Plots.scatter!(plt, [E1[1]], [E1[2]],
                   markersize=8, color=:blue, label="E₁")

    Plots.scatter!(plt, [E2[1]], [E2[2]],
                   markersize=8, color=:black, label="E₂")

    Plots.scatter!(plt, [E3[1]], [E3[2]],
                   markersize=8, color=:red, label="E₃")

    Plots.annotate!(plt, 500, 4000, Plots.text("E₁", 20, :blue))
    Plots.annotate!(plt, 2000, 700, Plots.text("E₂", 20, :black))
    Plots.annotate!(plt, 5500, 400, Plots.text("E₃", 20, :red))

    Plots.display(plt)
    Plots.savefig(plt, "phase_plane_p=2.png")

    return plt
end

run_phase_diagram_with_arrows()