# Julia code for the BLACK-AND-WHITE phenotype map highlighting the
# (S, S) bistable phenotype, for the SS05 parameter set of the model
# in Boareto et al. 2015 PNAS.
# Color-blind-safe rendering using patterns instead of color:
#
#   S       -> diagonal "/" hatch lines (continuous, broken at boundaries)
#   R       -> small filled dots
#   {S,R}   -> empty (pure white)       (unlikely to appear in this view)
#   {R,R}   -> solid black fill         (unlikely to appear in this view)
#   {S,S}   -> X-hatch (both "/" and "\" diagonals)
#
# Yuxin Liang, Ling Xue, Kun Zhao, Xiaoming Zheng

using CairoMakie
using LaTeXStrings
using Roots
_t_start = time()   # --- start runtime timer ---

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

# --- 3. Compute classification grid ---
res = 2000
d_min, d_max = 2e3, 1.5e4
log_n_min, log_n_max = 1.0, 5.0
d_vals = collect(range(d_min, d_max, length=res))
n_vals = collect((10.0).^range(log_n_min, log_n_max, length=res))
grid_data = zeros(Float64, length(n_vals), length(d_vals))

println("Calculating...")
for i in 1:length(n_vals), j in 1:length(d_vals)
    grid_data[i, j] = get_phenotype(d_vals[j], n_vals[i])
end

# --- 4. Base layer: white everywhere, black for {R,R} (rare here) ---
base_grid = ones(Float64, size(grid_data))
for k in eachindex(grid_data)
    if grid_data[k] == 3.0      # {R,R}
        base_grid[k] = 0.0
    end
end

# --- 5. Plotting ---
println("Plotting...")
begin
    fig = Figure(size = (1000, 800), fontsize = 30,
                 figure_padding = (15, 80, 15, 15))
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
              yscale = log10,
              aspect = 4/3)

    # Base layer: black for {R,R}, white elsewhere
    CairoMakie.heatmap!(ax, d_vals, n_vals, base_grid',
        colormap = [:black, :white], colorrange = (0, 1),
        interpolate = false)

    # ---- Helpers: y axis is log10, x axis is linear ----
    in_plot(d, log_n) =
        d_min <= d <= d_max && log_n_min <= log_n <= log_n_max
    function region_at(d, log_n)
        if !in_plot(d, log_n); return -1.0; end
        d_idx = clamp(round(Int,
            (d - d_min) / (d_max - d_min) * (res - 1) + 1), 1, res)
        n_idx = clamp(round(Int,
            (log_n - log_n_min) / (log_n_max - log_n_min) * (res - 1) + 1), 1, res)
        return grid_data[n_idx, d_idx]
    end

    # ---- Hatch lines (one direction or two) -----------------------------
    # Diagonals are drawn in the (D, log N) coordinate system.  In that
    # system, D spans ~1.3e4 and log10(N) spans 4 decades; the apparent
    # slope of the hatches is determined by the slope variable below.
    function hatch_in_region!(target_value, slope, spacing, linewidth)
        # Walk a family of lines  log_n = slope * (d - d_min) / d_step * log_step + c.
        # Compute the range of c that produces lines intersecting the plot.
        d_step = d_max - d_min
        log_step = log_n_max - log_n_min
        # At d=d_min, log_n = c.  At d=d_max, log_n = slope*log_step + c.
        # For the line to enter the plot we need c (or slope*log_step + c)
        # to lie within [log_n_min, log_n_max].
        c_min = log_n_min - max(0.0, slope) * log_step
        c_max = log_n_max - min(0.0, slope) * log_step
        c_step = spacing
        walk_d_step = (d_max - d_min) / 400.0
        for c in c_min:c_step:c_max
            seg_d = Float64[]
            seg_log_n = Float64[]
            d = d_min
            while d <= d_max
                # log_n is a linear function of (d - d_min)/(d_max - d_min)
                t = (d - d_min) / d_step
                log_n = slope * t * log_step + c
                if in_plot(d, log_n) && region_at(d, log_n) == target_value
                    push!(seg_d, d)
                    push!(seg_log_n, log_n)
                else
                    if length(seg_d) >= 2
                        lines!(ax, seg_d, 10.0 .^ seg_log_n;
                               color = :black, linewidth = linewidth)
                    end
                    empty!(seg_d); empty!(seg_log_n)
                end
                d += walk_d_step
            end
            if length(seg_d) >= 2
                lines!(ax, seg_d, 10.0 .^ seg_log_n;
                       color = :black, linewidth = linewidth)
            end
        end
    end

    # S region: forward diagonals only
    hatch_in_region!(0.0, +1.0, 0.30, 1.2)

    # {S,S} region: cross-hatch  (forward AND backward diagonals)
    hatch_in_region!(4.0, +1.0, 0.30, 1.2)
    hatch_in_region!(4.0, -1.0, 0.30, 1.2)

    # ---- R region: small filled dots (and {S,S} as "X" only if rare) ----
    pattern_res_d = 60
    pattern_res_n = 30
    pattern_d = collect(range(d_min, d_max, length=pattern_res_d))
    pattern_log_n = collect(range(log_n_min, log_n_max, length=pattern_res_n))
    pattern_n = 10.0 .^ pattern_log_n

    r_x, r_y = Float64[], Float64[]
    for di in 1:length(pattern_d), ni in 1:length(pattern_log_n)
        v = region_at(pattern_d[di], pattern_log_n[ni])
        if v == 1.0
            push!(r_x, pattern_d[di]); push!(r_y, pattern_n[ni])
        end
    end
    if !isempty(r_x)
        scatter!(ax, r_x, r_y; marker = :circle, markersize = 6, color = :black)
    end

    # ---- Region boundaries ----
    CairoMakie.contour!(ax, d_vals, n_vals, grid_data';
        levels    = [0.5, 1.5, 2.5, 3.5],
        color     = :black,
        linewidth = 2.0)

    # ---- Legend ----
    border_white = PolyElement(polycolor = :white,
                               polystrokecolor = :black, polystrokewidth = 2)
    border_black = PolyElement(polycolor = :black,
                               polystrokecolor = :black, polystrokewidth = 2)

    # Two parallel "/" lines reaching the boundary of the square
    s_line_1 = LineElement(points = Point2f[(0.0, 0.30), (0.70, 1.00)],
                           color = :black, linewidth = 1.5)
    s_line_2 = LineElement(points = Point2f[(0.30, 0.0), (1.00, 0.70)],
                           color = :black, linewidth = 1.5)

    # {S,S}: cross-hatch — the same two "/" lines, plus two "\" lines
    ss_line_1 = LineElement(points = Point2f[(0.0, 0.30), (0.70, 1.00)],
                            color = :black, linewidth = 1.5)
    ss_line_2 = LineElement(points = Point2f[(0.30, 0.0), (1.00, 0.70)],
                            color = :black, linewidth = 1.5)
    ss_line_3 = LineElement(points = Point2f[(0.0, 0.70), (0.70, 0.0)],
                            color = :black, linewidth = 1.5)
    ss_line_4 = LineElement(points = Point2f[(0.30, 1.00), (1.00, 0.30)],
                            color = :black, linewidth = 1.5)

    s_swatch  = [border_white, s_line_1, s_line_2]
    r_swatch  = [border_white,
                 MarkerElement(marker = :circle, markersize = 6, color = :black,
                               points = Point2f[(0.30, 0.30), (0.70, 0.30),
                                                (0.30, 0.70), (0.70, 0.70)])]
    ss_swatch = [border_white, ss_line_1, ss_line_2, ss_line_3, ss_line_4]

    Legend(fig,
        [s_swatch, r_swatch, ss_swatch],
        ["S", "R", "{S, S}"],
        "Phenotype";
        labelsize     = 38,
        titlesize     = 42,
        patchsize     = (55, 35),
        bbox          = ax.scene.viewport,
        halign        = :center,
        valign        = :top,
        width         = Relative(1.0),
        margin        = (0, 0, 10, 10),
        orientation   = :horizontal,
        nbanks        = 1,
        titleposition = :left,
        backgroundcolor = :white,
        framevisible  = true,
        framewidth    = 2.0)

    display(fig)
    save("phenotypebifur_SS05_BW.png", fig)
    println("Saved as phenotypebifur_SS05.png")
end

println("Running time: ", round(time() - _t_start, digits=2), " seconds")
