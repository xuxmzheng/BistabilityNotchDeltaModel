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
const Next_val = 12000
const I0       = 200.0

const N0, D0 = βN / γ, βD / γ
const K_const = (βN * kT) / (I0 * γI * γ)
const αd, αn = (kC * βD) / (γ^2), (kC * βN) / (γ^2)

# --- 2. Data Acquisition ---
d_vals = 10 .^ range(2, 10, length=800)

b_low_N, b_low_D = Point2f[], Point2f[]
b_mid_N, b_mid_D = Point2f[], Point2f[]
b_high_N, b_high_D = Point2f[], Point2f[]

for Dx in d_vals
    bi, bd, bn = K_const * Dx, (kT * Dx / γ) + 1.0, (kT * Next_val / γ) + 1.0
    f(N) = (2.0 - (1.0/(1.0 + (bi*N)^p))) - (αd*N*( (1.0/(1.0 + (bi*N)^p)) / (αn*N + bn) ) + bd*N)

    roots = Float64[]
    # Extremely dense test range to close the physical gap at the bifurcation point
    test_range = range(1e-12, 5.0, length=4000) 
    for i in 1:(length(test_range)-1)
        if f(test_range[i]) * f(test_range[i+1]) < 0
            push!(roots, find_zero(f, (test_range[i], test_range[i+1])))
        end
    end
    sort!(roots)

    if length(roots) == 1
        H = 1.0 / (1.0 + (bi * roots[1])^p)
        dr = H / (αn * roots[1] + bn)
        if Dx < 1e5 && roots[1] > 0.4
            push!(b_high_N, Point2f(Dx, roots[1] * N0)); push!(b_high_D, Point2f(Dx, dr * D0))
        else
            push!(b_low_N, Point2f(Dx, roots[1] * N0)); push!(b_low_D, Point2f(Dx, dr * D0))
        end
    elseif length(roots) >= 3
        for (i, nr) in enumerate(roots[1:3])
            H = 1.0 / (1.0 + (bi * nr)^p)
            dr = H / (αn * nr + bn)
            if i == 1
                push!(b_low_N, Point2f(Dx, nr * N0)); push!(b_low_D, Point2f(Dx, dr * D0))
            elseif i == 2
                push!(b_mid_N, Point2f(Dx, nr * N0)); push!(b_mid_D, Point2f(Dx, dr * D0))
            elseif i == 3
                push!(b_high_N, Point2f(Dx, nr * N0)); push!(b_high_D, Point2f(Dx, dr * D0))
            end
        end
    end
end

# --- 3. NaN Injection ---
function clean_branch(pts)
    if isempty(pts) return pts end
    cleaned = [pts[1]]
    for i in 2:length(pts)
        if abs(log10(pts[i][1]) - log10(pts[i-1][1])) > 0.02 
            push!(cleaned, Point2f(NaN, NaN))
        end
        push!(cleaned, pts[i])
    end
    return cleaned
end

b_low_N_c, b_low_D_c = clean_branch(b_low_N), clean_branch(b_low_D)
b_mid_N_c, b_mid_D_c = clean_branch(b_mid_N), clean_branch(b_mid_D)
b_high_N_c, b_high_D_c = clean_branch(b_high_N), clean_branch(b_high_D)

# --- 4. Plotting (Log-Log Scale) ---
begin
    fig = Figure(size = (1000, 800), fontsize = 30, figure_padding = (20, 85, 20, 20))
    
    # Added yscale = log10 here
    ax = Axis(fig[1, 1], 
              title = "Bifurcation Trace (p=$p, Next=$Next_val)", 
              xlabel = "Dext", ylabel = "molecules",
              xscale = log10, yscale = log10, aspect = 4/3,
              titlesize = 52,
              xlabelsize = 52, # Larger xlabel
              ylabelsize = 52, # Larger ylabel
              xticklabelsize = 52, 
              yticklabelsize = 52)           

    ax.xticks = ([10^2, 10^4, 10^6, 10^8, 10^10], ["10²", "10⁴", "10⁶", "10⁸", "10¹⁰"])
    ax.yticks = ([10^-2, 1, 10^2, 10^4, 10^6], ["10⁻²", "10⁰", "10²", "10⁴", "10⁶"])
    
    CairoMakie.xlims!(ax, 10^2, 10^10)
    CairoMakie.ylims!(ax, 10^-3, 10^7)

    lines!(ax, b_low_N_c,  color = :blue, linewidth = 4)
    lines!(ax, b_low_D_c,  color = :red,  linewidth = 4)
    lines!(ax, b_high_N_c, color = :blue, linewidth = 4)
    lines!(ax, b_high_D_c, color = :red,  linewidth = 4)
    
    lines!(ax, b_mid_N_c, color = (:blue, 0.6), linestyle = :dash, linewidth = 3)
    lines!(ax, b_mid_D_c, color = (:red, 0.6),  linestyle = :dash, linewidth = 3)

    leg_elements = [
        LineElement(color = :blue, linewidth = 4),
        LineElement(color = :red, linewidth = 4),
        LineElement(color = :black, linestyle = :dash, linewidth = 3)
    ]
    
Legend(fig, leg_elements, ["Notch", "Delta", "Unstable"], "",
       labelsize = 52, 
       titlesize = 52,
       patchsize = (50, 20), # Increase 100 to make the lines even longer
       bbox = ax.scene.viewport, 
       halign = :right, 
       valign = :top,            
       margin = (25, 25, 25, 25), 
       orientation = :horizontal,
       nbanks = 1, 
       titleposition = :left,      
       backgroundcolor = (:white, 0.85), 
       framevisible = true)

    resize_to_layout!(fig)
    display(fig)
    save("bifurcation_loglog_p=2_Next=12000.png", fig)
end