# Julia code for generating a BLACK-AND-WHITE zoom-in of the phenotype
# diagram near the quadruple point (p = 2) for the model in Boareto et al
# 2015 PNAS.
# Color-blind-safe rendering using patterns instead of color:
#
#   S       -> diagonal "/" hatch lines (continuous, broken at boundaries)
#   R       -> small filled dots
#   {S,R}   -> empty (pure white)
#   {R,R}   -> solid black fill
#   {S,S}   -> "X" marks  (very rare in this zoom; safety net)
#
# Yuxin Liang, Ling Xue, Kun Zhao, Xiaoming Zheng

using CairoMakie
using LaTeXStrings
using Roots
_t_start = time()   # --- start runtime timer ---

# --- 1. Parameters ---
const p = 2
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
    f(N) = (2.0 - (1.0/(1.0 + (bi*N)^p))) -
           (αd*N*( (1.0/(1.0 + (bi*N)^p)) / (αn*N + bn) ) + bd*N)
    roots = Float64[]
    test_range = range(1e-9, 2.1, length=300)
    for i in 1:(length(test_range)-1)
        try
            if f(test_range[i]) * f(test_range[i+1]) < 0
                push!(roots, find_zero(f, (test_range[i], test_range[i+1])))
            end
        catch; continue; end
    end
    stable_roots = length(roots) >= 3 ? [roots[1], roots[end]] :
                   (length(roots) > 0 ? [roots[1]] : Float64[])
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

# --- 3. Compute classification grid (linear zoom near the quadruple point) ---
res = 60
d_min, d_max = 661.5, 662.5
n_min, n_max = 8812.0, 8813.0
d_vals = collect(range(d_min, d_max, length=res))
n_vals = collect(range(n_min, n_max, length=res))
grid_data = zeros(Float64, length(n_vals), length(d_vals))

println("Calculating...")
for i in 1:length(n_vals), j in 1:length(d_vals)
    grid_data[i, j] = get_phenotype(d_vals[j], n_vals[i])
end

# --- 4. Build the base layer: white everywhere, black for {R,R} ---
base_grid = ones(Float64, size(grid_data))
for k in eachindex(grid_data)
    if grid_data[k] == 3.0          # {R,R}
        base_grid[k] = 0.0
    end
end

# --- 5. Plotting ---
println("Plotting...")
begin
    fig = Figure(size = (1000, 800), fontsize = 30, figure_padding = (15, 80, 15, 15))
    gl = fig[1, 1] = GridLayout()

    ax = Axis(gl[1, 1],
              title = L"\text{Cell Phenotypes } (p=%$p)",
              titlesize = 45,
              xlabel = L"\text{External Delta } (D_{ext})",
              ylabel = L"\text{External Notch } (N_{ext})",
              xlabelsize = 62,
              ylabelsize = 62,
              xticklabelsize = 50,
              yticklabelsize = 50,
              xscale = identity,
              yscale = identity,
              aspect = 4/3)

    CairoMakie.xlims!(ax, d_min, d_max)
    CairoMakie.ylims!(ax, n_min, n_max)

    # Base layer: black for {R,R}, white for everything else
    CairoMakie.heatmap!(ax, d_vals, n_vals, base_grid',
        colormap = [:black, :white], colorrange = (0, 1),
        interpolate = false)

    # ---- Helper: region lookup in linear coords ----
    in_plot(d, n) = d_min <= d <= d_max && n_min <= n <= n_max
    function region_at(d, n)
        if !in_plot(d, n); return -1.0; end
        d_idx = clamp(round(Int,
            (d - d_min) / (d_max - d_min) * (res - 1) + 1), 1, res)
        n_idx = clamp(round(Int,
            (n - n_min) / (n_max - n_min) * (res - 1) + 1), 1, res)
        return grid_data[n_idx, d_idx]
    end

    # ---- Continuous diagonal hatch lines in the S region ----
    # In the linear (D, N) zoom (each axis spans 1 unit) we use
    # N = D + c with slope 1, evenly spaced and clipped to the S region.
    hatch_spacing = 0.08      # spacing between hatch lines, in plot units
    walk_step     = 0.005     # sampling step along each line, in plot units
    c_min = n_min - d_max
    c_max = n_max - d_min

    for c in c_min:hatch_spacing:c_max
        seg_d = Float64[]
        seg_n = Float64[]
        d = d_min
        while d <= d_max
            n = d + c
            if in_plot(d, n) && region_at(d, n) == 0.0
                push!(seg_d, d)
                push!(seg_n, n)
            else
                if length(seg_d) >= 2
                    lines!(ax, seg_d, seg_n; color = :black, linewidth = 1.2)
                end
                empty!(seg_d); empty!(seg_n)
            end
            d += walk_step
        end
        if length(seg_d) >= 2
            lines!(ax, seg_d, seg_n; color = :black, linewidth = 1.2)
        end
    end

    # ---- R and {S,S} regions: marker overlays at a coarse grid ----
    pattern_res_d = 60
    pattern_res_n = 60
    pattern_d = collect(range(d_min, d_max, length=pattern_res_d))
    pattern_n = collect(range(n_min, n_max, length=pattern_res_n))

    r_x, r_y   = Float64[], Float64[]
    ss_x, ss_y = Float64[], Float64[]

    for di in 1:length(pattern_d), ni in 1:length(pattern_n)
        v = region_at(pattern_d[di], pattern_n[ni])
        if     v == 1.0;  push!(r_x,  pattern_d[di]); push!(r_y,  pattern_n[ni])
        elseif v == 4.0;  push!(ss_x, pattern_d[di]); push!(ss_y, pattern_n[ni])
        end
    end

    if !isempty(r_x)
        scatter!(ax, r_x, r_y; marker = :circle, markersize = 6, color = :black)
    end
    if !isempty(ss_x)
        scatter!(ax, ss_x, ss_y; marker = 'X', markersize = 18, color = :black)
    end

    # ---- Region boundaries (solid black contour lines) ----
    CairoMakie.contour!(ax, d_vals, n_vals, grid_data';
        levels    = [0.5, 1.5, 2.5, 3.5],
        color     = :black,
        linewidth = 2.0)

    # ---- Legend (sidebar to the right of the plot) ----
    border_white = PolyElement(polycolor = :white,
                               polystrokecolor = :black, polystrokewidth = 2)
    border_black = PolyElement(polycolor = :black,
                               polystrokecolor = :black, polystrokewidth = 2)

    s_line_1 = LineElement(
                   points = Point2f[(0.0, 0.30), (0.70, 1.00)],
                   color  = :black, linewidth = 1.5)
    s_line_2 = LineElement(
                   points = Point2f[(0.30, 0.0), (1.00, 0.70)],
                   color  = :black, linewidth = 1.5)

    s_swatch  = [border_white, s_line_1, s_line_2]
    r_swatch  = [border_white,
                 MarkerElement(marker = :circle, markersize = 6, color = :black,
                               points = Point2f[(0.30, 0.30), (0.70, 0.30),
                                                (0.30, 0.70), (0.70, 0.70)])]
    sr_swatch = [border_white]
    rr_swatch = [border_black]
    
    Legend(fig,
        [s_swatch, r_swatch, sr_swatch, rr_swatch],
        ["S", "R", "{S,R}", "{R,R}"],
        "Phenotype";

        bbox          = ax.scene.viewport,
        halign        = :center,
        valign        = :top,

        width         = Relative(1.0),
        tellwidth     = false,

        orientation   = :horizontal,
        nbanks        = 1,
        titleposition = :left,

        patchsize     = (32, 22),
        labelsize     = 26,
        titlesize     = 30,
        colgap        = 10,
        margin        = (4, 4, 4, 4),

        framevisible  = true,
        framewidth    = 1.5,
        backgroundcolor = :white
    )

    display(fig)
    save("pheno_p2_special_BW.png", fig)
    println("Saved as pheno_p2_special_BW.png")
end

println("Running time: ", round(time() - _t_start, digits=2), " seconds")
