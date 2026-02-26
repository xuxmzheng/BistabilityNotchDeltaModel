using CairoMakie
using Roots

# --- 1. Parameters ---
const p        = 2
const γ        = 0.1
const γI       = 0.5
const kC       = 5e-4
const kT       = 5e-5
const βN       = 500.0
const βD       = 1000.0
const Next_val = 8700
const I0       = 200.0

const N0, D0 = βN / γ, βD / γ
const K_const = (βN * kT) / (I0 * γI * γ)
const αd, αn = (kC * βD) / (γ^2), (kC * βN) / (γ^2)

# --- 2. Data Acquisition ---
d_vals = range(600, 700, length=800)

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
        log_y_jump = abs(log10(pts[i][2]) - log10(pts[i-1][2]))
        if abs(pts[i][1] - pts[i-1][1]) > 5.0 || log_y_jump > 0.5
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
    fig = Figure(size = (1000, 750), fontsize = 26, figure_padding = (20, 80, 20, 20))
    ax = Axis(fig[1, 1], 
              title = "Bifurcation Regions (p=$p, Next=$Next_val)", 
              xlabel = "Dext", ylabel = "molecules",
              yscale = log10,
              titlesize = 42, xlabelsize = 42, ylabelsize = 42,
              xticklabelsize = 32, yticklabelsize = 32,
              # --- REMOVE GRID LINES ---
              xgridvisible = false, 
              ygridvisible = false,
              xminorgridvisible = false,
              yminorgridvisible = false)  

    ax.yticks = ([1e-1, 1e0, 1e1, 1e2, 1e3, 1e4], ["10⁻¹", "10⁰", "10¹", "10²", "10³", "10⁴"])
    CairoMakie.ylims!(ax, 80, 12000) 
    CairoMakie.xlims!(ax, 600, 700)

    # Shading vertical regions
    CairoMakie.vspan!(ax, 600, 663, color = (:red, 0.1))    
    CairoMakie.vspan!(ax, 663, 676.5, color = (:yellow, 0.1)) 
    CairoMakie.vspan!(ax, 676.5, 700, color = (:lime, 0.1))   

    # Labels at the bottom
    text_y = 1200
    text!(ax, 625, text_y, text = "S", fontsize = 40, align = (:center, :center), color = :black)
    text!(ax, 670, text_y, text = "SR", fontsize = 40, align = (:center, :center), color = :black)
    text!(ax, 685.0, text_y, text = "RR", fontsize = 40, align = (:center, :center), color = :black)

    # Plot Branches
    lines!(ax, b1_N_c, color = :blue, linewidth = 4)
    lines!(ax, b1_D_c, color = :red,  linewidth = 4)
    lines!(ax, b2_N_c, color = :blue, linewidth = 2)
    lines!(ax, b2_D_c, color = :red,  linewidth = 2)
    lines!(ax, bm_N_c, color = (:blue, 0.4), linestyle = :dash, linewidth = 3)
    lines!(ax, bm_D_c, color = (:red, 0.4),  linestyle = :dash, linewidth = 3)

    # --- ADDED POINT ---
    # marker = :circle, markersize = 20, color = :black
    CairoMakie.scatter!(ax, [676.5], [920], color = :black, markersize = 50, marker = :circle)

    leg_elements = [LineElement(color = :blue, linewidth = 4),
                    LineElement(color = :red, linewidth = 4),
                    LineElement(color = :black, linestyle = :dash, linewidth = 3)]
    
    Legend(fig[1, 2], leg_elements, ["Notch", "Delta", "Unstable"],
           labelsize = 32, 
           patchsize = (40, 15),
           tellwidth = true,   
           tellheight = false) 

    display(fig)
    save("bifurcation_special_p2_lower.png", fig)
end