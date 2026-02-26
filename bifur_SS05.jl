using CairoMakie
using Roots

function run_bifurcation()

    γ, γI = 0.148, 0.456
    p = 2.04
    kC = 4.31e-04
    kT = 9.95e-05
    βN = 252.46
    βD = 1819.26
    I0 = 273.85
    Next_val = 100

    N0, D0 = βN / γ, βD / γ
    K_const = (βN * kT) / (I0 * γI * γ)
    αd, αn = (kC * βD) / (γ^2), (kC * βN) / (γ^2)

    d_vals = range(5000.0, 15000.0, length=800)

    b_low_N,  b_low_D  = Point2f[], Point2f[]
    b_mid_N,  b_mid_D  = Point2f[], Point2f[]
    b_high_N, b_high_D = Point2f[], Point2f[]

    prev_root_count = 0

    for Dx in d_vals

        bi = K_const * Dx
        bd = (kT * Dx / γ) + 1.0
        bn = (kT * Next_val / γ) + 1.0

        f(N) = (2.0 - (1.0/(1.0 + (bi*N)^p))) -
               (αd*N*( (1.0/(1.0 + (bi*N)^p)) / (αn*N + bn) ) + bd*N)

        roots = Float64[]
        test_range = range(1e-12, 15.0, length=3000)

        for i in 1:length(test_range)-1
            if f(test_range[i]) * f(test_range[i+1]) < 0
                push!(roots, find_zero(f, (test_range[i], test_range[i+1])))
            end
        end

        sort!(roots)

        # ---- Deduplicate (fold fix) ----
        unique_roots = Float64[]
        tol = 1e-4
        for r in roots
            if isempty(unique_roots) || abs(r - unique_roots[end]) > tol
                push!(unique_roots, r)
            end
        end
        roots = unique_roots

        isempty(roots) && continue

        # ---- If root count changes, break curves ----
        if prev_root_count != length(roots)
            push!(b_low_N,  Point2f(NaN, NaN))
            push!(b_low_D,  Point2f(NaN, NaN))
            push!(b_mid_N,  Point2f(NaN, NaN))
            push!(b_mid_D,  Point2f(NaN, NaN))
            push!(b_high_N, Point2f(NaN, NaN))
            push!(b_high_D, Point2f(NaN, NaN))
        end

        # ---- Monostable ----
        if length(roots) == 1
            nr = roots[1]
            H = 1.0 / (1.0 + (bi * nr)^p)
            dr = H / (αn * nr + bn)

            push!(b_low_N, Point2f(Dx, nr * N0))
            push!(b_low_D, Point2f(Dx, dr * D0))

        # ---- Bistable ----
        elseif length(roots) >= 3

            # Low
            nr = roots[1]
            H = 1.0 / (1.0 + (bi * nr)^p)
            dr = H / (αn * nr + bn)
            push!(b_low_N, Point2f(Dx, nr * N0))
            push!(b_low_D, Point2f(Dx, dr * D0))

            # Mid (unstable)
            nr = roots[2]
            H = 1.0 / (1.0 + (bi * nr)^p)
            dr = H / (αn * nr + bn)
            push!(b_mid_N, Point2f(Dx, nr * N0))
            push!(b_mid_D, Point2f(Dx, dr * D0))

            # High
            nr = roots[3]
            H = 1.0 / (1.0 + (bi * nr)^p)
            dr = H / (αn * nr + bn)
            push!(b_high_N, Point2f(Dx, nr * N0))
            push!(b_high_D, Point2f(Dx, dr * D0))
        end

        prev_root_count = length(roots)
    end

    # ---- Plot ----
    fig = Figure(size=(1000,750))
    # ax = Axis(fig[1,1], xlabel="Dext", ylabel="molecules")

    ax = Axis(fig[1, 1], 
        title = "Bifurcation (p=$p, Next=$Next_val)", 
        xlabel = "Dext", ylabel = "molecules",
        titlesize = 52,
        xlabelsize = 52, # Larger xlabel
        ylabelsize = 52, # Larger ylabel
        xticklabelsize = 32, # Larger tick numbers
        yticklabelsize = 42) # Larger tick numbers

    lines!(ax, b_low_N,  color=:blue, linewidth=4)
    lines!(ax, b_high_N, color=:blue, linewidth=4)
    lines!(ax, b_low_D,  color=:red,  linewidth=4)
    lines!(ax, b_high_D, color=:red,  linewidth=4)

    lines!(ax, b_mid_N, color=(:blue,0.4), linestyle=:dash, linewidth=3)
    lines!(ax, b_mid_D, color=(:red,0.4),  linestyle=:dash, linewidth=3)

    display(fig)
    save("bifur_SS05.png", fig)
end

run_bifurcation()