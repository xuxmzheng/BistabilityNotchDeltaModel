using Random
using ProgressMeter

# ----------------------------
# 1. Global Constants
# ----------------------------
const N_TRIALS = Int(1e7)         
const P_MODE = "random"           
const P_FIXED = 1.5
const P_MIN, P_MAX = 1.05, 1.95

const LOG10_N_MIN = -8.0
const LOG10_N_MAX = 8.0
const GRID_PER_DECADE = 80      

const BISECTION_ITERS = 80

const ALPHA_N_RANGE = (1e-5, 100.0)
const ALPHA_D_RANGE = (1e-5, 100.0)  
const BETA_I_RANGE  = (1e-5, 2500.0)

# ----------------------------
# 2. Optimized Root Finding
# ----------------------------
function bisection_root(f, lo, hi, iters=BISECTION_ITERS)
    flo, fhi = f(lo), f(hi)
    if !isfinite(flo) || !isfinite(fhi) || flo * fhi > 0 return nothing end
    for _ in 1:iters
        mid = 0.5 * (lo + hi)
        fmid = f(mid)
        if !isfinite(fmid) return nothing end
        if flo * fmid <= 0
            hi, fhi = mid, fmid
        else
            lo, flo = mid, fmid
        end
    end
    return 0.5 * (lo + hi)
end

function count_positive_roots_fast(p, alpha_n, alpha_d, beta_n, beta_d, beta_i, Ngrid)
    c_p1 = (beta_d * beta_n - 2.0 * alpha_n) / (alpha_n * beta_d)
    c_p0 = -(2.0 * beta_n) / (alpha_n * beta_d)
    inv_beta_ip = 1.0 / (beta_i ^ p)
    c_2, c_1, c_0 = inv_beta_ip, (alpha_d + beta_n * beta_d - alpha_n) / (alpha_n * beta_d) * inv_beta_ip, -(beta_n) / (alpha_n * beta_d) * inv_beta_ip

    Np = Ngrid .^ p
    N2 = Ngrid .* Ngrid
    vals = Np .* (N2 .+ c_p1 .* Ngrid .+ c_p0) .+ (c_2 .* N2 .+ c_1 .* Ngrid .+ c_0)

    s = sign.(vals)
    change_indices = findall(x -> x != 0, diff(s))
    
    roots = Float64[]
    for idx in change_indices
        a, b = Ngrid[idx], Ngrid[idx+1]
        f_loc(N) = (N^p * (N^2 + c_p1*N + c_p0) + (c_2*N^2 + c_1*N + c_0))
        r = bisection_root(f_loc, a, b)
        if !isnothing(r) push!(roots, r) end
    end
    return roots
end

# ----------------------------
# 3. Main Loop
# ----------------------------
function run_trials(; n_trials=N_TRIALS)
    Random.seed!(0)
    rng = Random.MersenneTwister(12345)
    best = 0
    best_cases = []
    
    Ngrid = 10.0 .^ range(LOG10_N_MIN, LOG10_N_MAX, length=Int(round((LOG10_N_MAX-LOG10_N_MIN)*GRID_PER_DECADE))+1)
    
    # Cyan meter with [=> ] glyphs
    p_bar = Progress(n_trials; dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=40, color=:cyan)
    
    for t in 1:n_trials
        param_p = P_MODE == "fixed" ? P_FIXED : rand(rng, P_MIN:P_MAX)
        
        # Parameter Sampling
        an, ad = rand(rng, ALPHA_N_RANGE[1]:ALPHA_N_RANGE[2]), rand(rng, ALPHA_D_RANGE[1]:ALPHA_D_RANGE[2])
        L, U = an - ad, 2.0 * an
        if U > 0 && L < U
            prod = rand(rng, max(L, 1e-6):(0.999999 * U))
            bn = rand(rng, (1.0 + 1e-6):8.0) 
            bd = prod / bn
            if bd > 1.0
                bi = rand(rng, BETA_I_RANGE[1]:BETA_I_RANGE[2])
                
                # Fast Vectorized Core
                roots = count_positive_roots_fast(param_p, an, ad, bn, bd, bi, Ngrid)
                k = length(roots)

                if k > best
                    best = k
                    best_cases = [(param_p, an, ad, bn, bd, bi, roots)]
                    # Only report if 5 roots are found
                    if k >= 5
                        print("\r\u1b[K") # Clear the bar line
                        println("\n[Trial $t] SUCCESS: Found $k roots!")
                        println("Parameters: p=$param_p, alpha_n=$an, alpha_d=$ad, beta_n=$bn, beta_d=$bd, beta_i=$bi")
                    end
                elseif k == best && best > 0 && length(best_cases) < 5
                    push!(best_cases, (param_p, an, ad, bn, bd, bi, roots))
                end
            end
        end

        # Update bar. Semicolon ensures no return value leakage.
        next!(p_bar);

        if best >= 5
            finish!(p_bar)
            break
        end
    end
    
    return (best, best_cases)
end

# ----------------------------
# 4. Execution
# ----------------------------
# Trailing semicolon is mandatory to suppress the final Result metadata
res_best, res_cases = run_trials();

println("\n" * "="^40)
println("Search finished. Max roots found: $res_best")
if res_best > 0
    println("Parameters of first best case:")
    println(res_cases[1][1:6])
end
println("="^40)