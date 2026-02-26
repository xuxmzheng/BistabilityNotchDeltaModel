using CairoMakie
using Roots
using Plots

#    Params: bN=252.46, bD=1819.26, g=0.148, gI=0.456, kC=4.31e-04, kT=9.95e-05, 
#    I0=273.85, Dext=7185.64, Next=573.19, p=2.04

# --- 1. Parameters ---
const γ, γI = 0.148, 0.456
const p = 2.04
const kC = 4.31e-04
const kT = 9.95e-05
const βN = 252.46
const βD = 1819.26
const I0 = 273.85

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
res = 500 
d_vals = collect(range(2e3, 1.5e4, length=res)) 
n_vals = collect((10.0).^range(1, 5, length=res))
grid_data = zeros(Float64, length(n_vals), length(d_vals))

println("Calculating...")
for i in 1:length(n_vals), j in 1:length(d_vals)
    grid_data[i, j] = get_phenotype(d_vals[j], n_vals[i])
end

# --- 4. Plotting (ULTRA LARGE FONTS) ---
begin
    # Massive canvas size (2000x1400) and base fontsize 60
    fig = Figure(size = (2000, 1400), fontsize = 60, figure_padding = (40, 120, 40, 40))
    
    # alignmode=Outside() ensures text doesn't overlap the plot area
    ax = Axis(fig[1, 1], 
              title = "Cell Phenotypes", 
              titlesize = 80, 
              xlabel = "External Delta (Dext)", 
              ylabel = "External Notch (Next)",
              xlabelsize = 70,
              ylabelsize = 70,
              xticklabelsize = 55,
              yticklabelsize = 55,
              xticksize = 25, 
              yticksize = 25,
              xtickwidth = 4,
              ytickwidth = 4,
              yscale = log10,
              alignmode = Outside()) 

    cmap = [:red, :blue, :cyan]

    CairoMakie.heatmap!(ax, d_vals, n_vals, grid_data', 
              colormap = cmap, 
              colorrange = (0, 2), 
              interpolate = false)

    labels = ["S", "R", "{S, S}"]
    elements = [PolyElement(polycolor = cmap[i]) for i in 1:3]
    
    Legend(fig[1, 2], elements, labels, "Phenotype",
           labelsize = 55,
           titlesize = 65,
           titlegap = 30,
           padding = (30, 30, 30, 30),
           patchsize = (80, 80)) # Massive legend boxes

    display(fig)
    save("phenotypebifur_SS05.png", fig)
end