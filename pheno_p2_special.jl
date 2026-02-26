# Julia code for generating phenotype bifurcation w.r.t Dext and Next
# for the model in Boareto et al 2015 PNAS
# Yuxin Liang, Ling Xue, Kun Zhao, Xiaoming Zheng
# January 2026

using CairoMakie
using Roots

# --- 1. Parameters ---
const p = 2.0
const γ = 0.1
const γI = 0.5
const kC = 5e-4
const kT = 5e-5
const βN = 500.0
const βD = 1000.0
const I0 = 200.0

const N0, D0 = βN / γ, βD / γ
const K_const = (βN * kT) / (I0 * γI * γ)
const αd, αn = (kC * βD) / (γ^2), (kC * βN) / (γ^2)

# --- 2. Robust Classification ---
function get_phenotype(Dext::Float64, Next::Float64)
    bi = K_const * Dext
    bd = (kT * Dext / γ) + 1.0
    bn = (kT * Next / γ) + 1.0
    f(N) = (2.0 - (1.0/(1.0 + (bi*N)^p))) - (αd*N*( (1.0/(1.0 + (bi*N)^p)) / (αn*N + bn) ) + bd*N)
    roots = Float64[]
    test_range = range(1e-9, 2.1, length=300)
    for i in 1:(length(test_range)-1)
        try
            if f(test_range[i]) * f(test_range[i+1]) < 0
                push!(roots, find_zero(f, (test_range[i], test_range[i+1])))
            end
        catch; continue; end
    end
    stable_roots = length(roots) >= 3 ? [roots[1], roots[end]] : (length(roots) > 0 ? [roots[1]] : Float64[])
    if isempty(stable_roots); return 0.0; end
    states = Symbol[]
    for nr in stable_roots
        H = 1.0 / (1.0 + (bi * nr)^p)
        dr = H / (αn * nr + bn)
        push!(states, (nr * N0 > dr * D0) ? :R : :S)
    end
    if length(states) == 1
        return states[1] == :S ? 0.0 : 1.0
    else
        if all(s -> s == :R, states); return 3.0 
        elseif all(s -> s == :S, states); return 4.0 
        else; return 2.0; end
    end
end

# --- 3. Compute Grid ---
res = 600 
# Reverted to linear ranges
d_vals = collect(range(661.5, 662.5, length=res)) 
n_vals = collect(range(8812, 8813, length=res))
grid_data = zeros(Float64, length(n_vals), length(d_vals))

println("Calculating...")
for i in 1:length(n_vals), j in 1:length(d_vals)
    grid_data[i, j] = get_phenotype(d_vals[j], n_vals[i])
end

# --- 4. Plotting (Legend Outside to the Right) ---
begin
    fig = Figure(size = (1200, 800), fontsize = 30, figure_padding = 15)
    
    gl = fig[1, 1] = GridLayout()

    # Axis remains in the main grid cell
    ax = Axis(gl[1, 1], 
              title = "Cell Phenotypes (p=$p)", 
              titlesize = 50, 
              xlabel = "External Delta (Dext)", 
              ylabel = "External Notch (Next)",
              xlabelsize = 50,
              ylabelsize = 50,
              xticklabelsize = 40,
              yticklabelsize = 40,
              xscale = identity, 
              yscale = identity,
              aspect = 4/3) 
    
    CairoMakie.ylims!(ax, 8812, 8813) 
    CairoMakie.xlims!(ax, 661.5, 662.5)

    cmap = [:red, :blue, :yellow, :lime, :cyan]
    CairoMakie.heatmap!(ax, d_vals, n_vals, grid_data', colormap = cmap, colorrange = (0, 4), interpolate = false)

    labels = ["S", "R", "{S, R}", "{R, R}"]
    elements = [PolyElement(polycolor = cmap[i]) for i in 1:4]
    
    # Legend placed in the second column (gl[1, 2])
    # Orientation switched to :vertical for a better side-bar look
    Legend(gl[1, 2], elements, labels, "",
           labelsize = 28,           
           titlesize = 32,           
           patchsize = (40, 40),     
           halign = :left,                    
           orientation = :vertical,
           nbanks = 1,                
           titleposition = :top,     
           tellheight = false,
           tellwidth = true,         # Reserves horizontal space for the legend
           backgroundcolor = (:white, 0.75),
           framevisible = true,
           framewidth = 2.0)

    colgap!(gl, 20) # Space between the plot and the legend
    resize_to_layout!(fig)

    display(fig)
    save("pheno_p2_special.png", fig)
end