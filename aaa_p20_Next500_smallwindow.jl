using CairoMakie
using Roots

# --- 1. Parameters ---
const p        = 20
const γ        = 0.1
const γI       = 0.5
const kC       = 5e-4
const kT       = 5e-5
const βN       = 500.0
const βD       = 1000.0
const Next_val = 500
const I0       = 200.0

const N0, D0 = βN / γ, βD / γ
const K_const = (βN * kT) / (I0 * γI * γ)
const αd, αn = (kC * βD) / (γ^2), (kC * βN) / (γ^2)

# --- 2. Data Acquisition ---
# Changed 3000.0 to 10000.0 to match your plot limits
d_vals = range(0.0, 10000.0, length=800)

b1_N, b1_D = Point2f[], Point2f[] 
b2_N, b2_D = Point2f[], Point2f[] 
b_mid_N, b_mid_D = Point2f[], Point2f[] 

for Dx in d_vals
    bi, bd, bn = K_const * Dx, (kT * Dx / γ) + 1.0, (kT * Next_val / γ) + 1.0
    f(N) = (2.0 - (1.0/(1.0 + (bi*N)^p))) - (αd*N*( (1.0/(1.0 + (bi*N)^p)) / (αn*N + bn) ) + bd*N)

    roots = Float64[]
    test_range = range(1e-12, 15.0, length=2000) 
    for i in 1:(length(test_range)-1)
        if f(test_range[i]) * f(test_range[i+1]) < 0
            push!(roots, find_zero(f, (test_range[i], test_range[i+1])))
        end
    end
    sort!(roots)

    if length(roots) == 1
        nr = roots[1]
        H = 1.0 / (1.0 + (bi * nr)^p)
        dr = H / (αn * nr + bn)
        
        if nr < 0.5
            push!(b2_N, Point2f(Dx, nr * N0)); push!(b2_D, Point2f(Dx, dr * D0))
        else
            push!(b1_N, Point2f(Dx, nr * N0)); push!(b1_D, Point2f(Dx, dr * D0))
        end
        
    elseif length(roots) >= 3
        for (i, nr) in enumerate(roots[1:3])
            H = 1.0 / (1.0 + (bi * nr)^p)
            dr = H / (αn * nr + bn)
            if i == 1 
                push!(b2_N, Point2f(Dx, nr * N0)); push!(b2_D, Point2f(Dx, dr * D0))
            elseif i == 2 
                push!(b_mid_N, Point2f(Dx, nr * N0)); push!(b_mid_D, Point2f(Dx, dr * D0))
            elseif i == 3 
                push!(b1_N, Point2f(Dx, nr * N0)); push!(b1_D, Point2f(Dx, dr * D0))
            end
        end
    end
end

# --- 3. Cleaning ---
function clean_branch(pts)
    if isempty(pts) return pts end
    cleaned = [pts[1]]
    for i in 2:length(pts)
        # Increased 5.0 to 20.0 because the gap between points is now larger (~12.5)
        if abs(pts[i][1] - pts[i-1][1]) > 20.0 || abs(pts[i][2] - pts[i-1][2]) > 500.0
            push!(cleaned, Point2f(NaN, NaN))
        end
        push!(cleaned, pts[i])
    end
    return cleaned
end

b1_N_c, b1_D_c = clean_branch(b1_N), clean_branch(b1_D)
b2_N_c, b2_D_c = clean_branch(b2_N), clean_branch(b2_D)
bm_N_c, bm_D_c = clean_branch(b_mid_N), clean_branch(b_mid_D)

# --- 4. Plotting ---
begin
    fig = Figure(size = (1000, 750), fontsize = 26, figure_padding = (20, 50, 20, 20))
    ax = Axis(fig[1, 1], 
              title = "Bifurcation (p=$p, Next=$Next_val)", 
              xlabel = "Dext", ylabel = "molecules",
              titlesize = 52,
              xlabelsize = 52, 
              ylabelsize = 52, 
              xticklabelsize = 36, 
              yticklabelsize = 36) 

    CairoMakie.xlims!(ax, 0, 10000)
    CairoMakie.ylims!(ax, 0, 9000)

    lines!(ax, b1_N_c, color = :blue, linewidth = 4)
    lines!(ax, b1_D_c, color = :red,  linewidth = 4)
    lines!(ax, b2_N_c, color = :blue, linewidth = 4)
    lines!(ax, b2_D_c, color = :red,  linewidth = 4)
    lines!(ax, bm_N_c, color = (:blue, 0.4), linestyle = :dash, linewidth = 3)
    lines!(ax, bm_D_c, color = (:red, 0.4),  linestyle = :dash, linewidth = 3)

    leg_elements = [LineElement(color = :blue, linewidth = 4),
                    LineElement(color = :red, linewidth = 4),
                    LineElement(color = :black, linestyle = :dash, linewidth = 3)]
    
    Legend(fig[1, 1], leg_elements, ["Notch", "Delta", "Unstable"],
           labelsize = 52, 
           patchsize = (50, 20), 
           halign = :right, valign = :top,
           tellwidth = false, tellheight = false,
           margin = (10, 10, 10, 10),
           backgroundcolor = (:white, 0.8))

    display(fig)
    save("bifurcation_p20_Next500_smallwindow.png", fig)
end