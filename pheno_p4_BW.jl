# Julia code for generating a BLACK-AND-WHITE phenotype bifurcation
# w.r.t Dext and Next for the model in Boareto et al 2015 PNAS.
# Color-blind-safe rendering using patterns instead of color:
#
#   S       -> diagonal "/" marks
#   R       -> small filled dots
#   {S,R}   -> empty (pure white)
#   {R,R}   -> solid black fill
#   {S,S}   -> "X" marks  (rare for p=4, included as a safety net)
#
# Yuxin Liang, Ling Xue, Kun Zhao, Xiaoming Zheng

using CairoMakie
using LaTeXStrings
using Roots
_t_start = time()   # --- start runtime timer ---

# --- 1. Parameters ---
const p = 4
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

# --- 3. Compute high-resolution classification grid ---
res = 2000
d_vals = collect((10.0).^range(2, 10, length=res))
n_vals = collect((10.0).^range(2,  6, length=res))
grid_data = zeros(Float64, length(n_vals), length(d_vals))

println("Calculating high-resolution grid...")
for i in 1:length(n_vals), j in 1:length(d_vals)
    grid_data[i, j] = get_phenotype(d_vals[j], n_vals[i])
end

# --- 4. Build the base layer:
#         white everywhere except {R,R} which is filled black ---
base_grid = ones(Float64, size(grid_data))   # 1.0 = white
for k in eachindex(grid_data)
    if grid_data[k] == 3.0          # {R,R}
        base_grid[k] = 0.0          # 0.0 = black
    end
end

# --- 5. Plotting ---
println("Plotting...")
begin
    fig = Figure(size = (1000, 800), fontsize = 30,
                 figure_padding = (15, 40, 15, 15))
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
              xscale = log10,
              yscale = log10,
              aspect = 4/3)

    ax.xticks = ([10^2, 10^4, 10^6, 10^8, 10^10],
                 ["10²", "10⁴", "10⁶", "10⁸", "10¹⁰"])
    ax.yticks = ([10^2, 10^4, 10^6], ["10²", "10⁴", "10⁶"])
    CairoMakie.xlims!(ax, 10^2, 10^10)
    CairoMakie.ylims!(ax, 10^2, 10^6)

    # Base: black for {R,R}, white for everything else
    CairoMakie.heatmap!(ax, d_vals, n_vals, base_grid',
        colormap = [:black, :white], colorrange = (0, 1),
        interpolate = false)

    # ---- Helper: test whether a (log_d, log_n) point lies in a given region ----
    log_d_min, log_d_max = 2.0, 10.0
    log_n_min, log_n_max = 2.0,  6.0
    in_plot(log_d, log_n) =
        log_d_min <= log_d <= log_d_max && log_n_min <= log_n <= log_n_max
    function region_at(log_d, log_n)
        if !in_plot(log_d, log_n); return -1.0; end
        d_idx = clamp(round(Int,
            (log_d - log_d_min) / (log_d_max - log_d_min) * (res - 1) + 1), 1, res)
        n_idx = clamp(round(Int,
            (log_n - log_n_min) / (log_n_max - log_n_min) * (res - 1) + 1), 1, res)
        return grid_data[n_idx, d_idx]
    end

    # ---- Continuous diagonal hatch lines in the S region ----
    # In log-coordinates, the diagonal family is log(N) = log(D) + c.
    # Each line is sampled finely and split into segments wherever it
    # crosses out of the S region, so segments run from one region
    # boundary to the next.
    hatch_spacing = 0.55     # spacing between hatch lines, in log10 units
    walk_step     = 0.05     # sampling step along each line, in log10 units
    c_min = log_n_min - log_d_max
    c_max = log_n_max - log_d_min

    for c in c_min:hatch_spacing:c_max
        seg_log_d = Float64[]
        seg_log_n = Float64[]
        log_d = log_d_min
        while log_d <= log_d_max
            log_n = log_d + c
            if in_plot(log_d, log_n) && region_at(log_d, log_n) == 0.0
                push!(seg_log_d, log_d)
                push!(seg_log_n, log_n)
            else
                # End of an S segment: emit it
                if length(seg_log_d) >= 2
                    lines!(ax, 10.0 .^ seg_log_d, 10.0 .^ seg_log_n;
                           color = :black, linewidth = 1.2)
                end
                empty!(seg_log_d); empty!(seg_log_n)
            end
            log_d += walk_step
        end
        if length(seg_log_d) >= 2
            lines!(ax, 10.0 .^ seg_log_d, 10.0 .^ seg_log_n;
                   color = :black, linewidth = 1.2)
        end
    end

    # ---- R and {S,S} regions: marker overlays at a coarse grid ----
    pattern_res_d = 60
    pattern_res_n = 30
    pattern_d = collect((10.0).^range(2, 10, length=pattern_res_d))
    pattern_n = collect((10.0).^range(2,  6, length=pattern_res_n))

    r_x, r_y   = Float64[], Float64[]
    ss_x, ss_y = Float64[], Float64[]

    for di in 1:length(pattern_d), ni in 1:length(pattern_n)
        v = region_at(log10(pattern_d[di]), log10(pattern_n[ni]))
        if     v == 1.0;  push!(r_x,  pattern_d[di]); push!(r_y,  pattern_n[ni])
        elseif v == 4.0;  push!(ss_x, pattern_d[di]); push!(ss_y, pattern_n[ni])
        end
    end

    # R: small filled dots
    if !isempty(r_x)
        scatter!(ax, r_x, r_y; marker = :circle, markersize = 6, color = :black)
    end
    # {S,S}: "X" character (rare for p=4)
    if !isempty(ss_x)
        scatter!(ax, ss_x, ss_y; marker = 'X', markersize = 18, color = :black)
    end

    # ---- Region boundaries (solid black contour lines) ----
    # Phenotype values are 0 (S), 1 (R), 2 ({S,R}), 3 ({R,R}), 4 ({S,S});
    # half-integer levels separate the regions.
    CairoMakie.contour!(ax, d_vals, n_vals, grid_data';
        levels    = [0.5, 1.5, 2.5, 3.5],
        color     = :black,
        linewidth = 2.0)

    # ---- Legend ----
    # Each entry is a black-bordered square; the pattern is drawn inside.
    border_white = PolyElement(polycolor = :white,
                               polystrokecolor = :black, polystrokewidth = 2)
    border_black = PolyElement(polycolor = :black,
                               polystrokecolor = :black, polystrokewidth = 2)

    # Two parallel slanted lines that reach the boundary of the legend square.
    # Coordinates are in the swatch's local [0,1] x [0,1] frame.
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
        ["S", "R", "{S, R}", "{R, R}"],
        "Phenotype";
        labelsize     = 38,
        titlesize     = 42,
        patchsize     = (55, 35),
        bbox          = ax.scene.viewport,
        halign        = :center,
        valign        = :top,
        width         = Relative(0.99),
        margin        = (0, 0, 10, 10),
        orientation   = :horizontal,
        nbanks        = 1,
        titleposition = :left,
        backgroundcolor = :white,
        framevisible  = true,
        framewidth    = 2.0)

    resize_to_layout!(fig)
    display(fig)
    save("pheno_p4_BW.png", fig)
    println("Saved as pheno_p4_BW.png")
end

println("Running time: ", round(time() - _t_start, digits=2), " seconds")
