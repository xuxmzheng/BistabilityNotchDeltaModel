# Julia code: Parameter sweep with sorted descending phenotype report
# Hybrid Log-Uniform and Uniform sampling based on Boareto/Sprinzak biological ranges
# =====================================================================================
# GLOBAL ENVIRONMENT DISTRIBUTION (Across all 1000000 valid samples)
# Next > Dext: 499801     (49.98%)
# Next < Dext: 500199     (50.02%)
# -------------------------------------------------------------------------------------
# PHENOTYPE                      | TOTAL %    | Next > Dext %      | Next < Dext %
# -------------------------------------------------------------------------------------
# Monostable {R}                 | 66.76     % | 50.74             % | 49.26             %      
# Monostable {S}                 | 29.56     % | 49.72             % | 50.28             %      
# Bistable {R,S}                 | 3.34      % | 36.91             % | 63.09             %      
# Bistable {R,R}                 | 0.29      % | 60.77             % | 39.23             %      
# Bistable {S,S}                 | 0.05      % | 9.96              % | 90.04             %      
# =====================================================================================

using Roots, Printf, Random, ProgressMeter, UUIDs

# --- 1. Sweep Function ---
function run_robust_parameter_sweep(samples=100000)
    new_seed = hash((time_ns(), getpid(), uuid4()))
    Random.seed!(new_seed)
    
    println("Run started with unique seed: ", new_seed)

    results = Dict{String, Int}()
    env_counts = Dict{String, Dict{Symbol, Int}}() 
    
    total_next_high = 0
    total_dext_high = 0
    
    valid_count = 0
    ss_export_data = []

    log_ranges = Dict(
        :betaN  => (log10(50.0),   log10(10000.0)),
        :betaD  => (log10(50.0),   log10(10000.0)),
        :kC     => (log10(1e-6),   log10(1e-2)),
        :kT     => (log10(1e-6),   log10(1e-2)),
        :Dext   => (log10(50.0),   log10(10000.0)),
        :Next   => (log10(50.0),   log10(10000.0))
    )

    uni_ranges = Dict(
        :gamma  => (0.05, 0.2),
        :gammaI => (0.25, 1.0),
        :I0     => (100.0, 300.0),
        :p      => (1.0, 6.0)
    )

    p_bar = Progress(samples; dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=40, color=:cyan)

    for i in 1:samples
        c_betaN  = 10.0^(rand() * (log_ranges[:betaN][2] - log_ranges[:betaN][1]) + log_ranges[:betaN][1])
        c_betaD  = 10.0^(rand() * (log_ranges[:betaD][2] - log_ranges[:betaD][1]) + log_ranges[:betaD][1])
        c_kC     = 10.0^(rand() * (log_ranges[:kC][2] - log_ranges[:kC][1]) + log_ranges[:kC][1])
        c_kT     = 10.0^(rand() * (log_ranges[:kT][2] - log_ranges[:kT][1]) + log_ranges[:kT][1])
        c_Dext   = 10.0^(rand() * (log_ranges[:Dext][2] - log_ranges[:Dext][1]) + log_ranges[:Dext][1])
        c_Next   = 10.0^(rand() * (log_ranges[:Next][2] - log_ranges[:Next][1]) + log_ranges[:Next][1])
        
        c_gamma  = rand() * (uni_ranges[:gamma][2] - uni_ranges[:gamma][1]) + uni_ranges[:gamma][1]
        c_gammaI = rand() * (uni_ranges[:gammaI][2] - uni_ranges[:gammaI][1]) + uni_ranges[:gammaI][1]
        c_I0     = rand() * (uni_ranges[:I0][2] - uni_ranges[:I0][1]) + uni_ranges[:I0][1]
        c_p      = rand() * (uni_ranges[:p][2] - uni_ranges[:p][1]) + uni_ranges[:p][1]

        N0, D0 = c_betaN / c_gamma, c_betaD / c_gamma
        K_const = (c_betaN * c_kT) / (c_I0 * c_gammaI * c_gamma)
        αd, αn = (c_kC * c_betaD) / (c_gamma^2), (c_kC * c_betaN) / (c_gamma^2)
        
        bi = K_const * c_Dext
        bd = (c_kT * c_Dext / c_gamma) + 1.0
        bn = (c_kT * c_Next / c_gamma) + 1.0

        f(N) = (2.0 - (1.0/(1.0 + (bi*N)^c_p))) - (αd*N*( (1.0/(1.0 + (bi*N)^c_p)) / (αn*N + bn) ) + bd*N)
        
        roots = Float64[]
        test_range = range(1e-9, 2.1, length=300) 
        for j in 1:(length(test_range)-1)
            try
                if f(test_range[j]) * f(test_range[j+1]) < 0
                    push!(roots, find_zero(f, (test_range[j], test_range[j+1])))
                end
            catch; continue; end
        end

        stable_roots = length(roots) >= 3 ? [roots[1], roots[end]] : (length(roots) > 0 ? [roots[1]] : Float64[])
        
        if !isempty(stable_roots)
            valid_count += 1
            if c_Next > c_Dext
                total_next_high += 1
            else
                total_dext_high += 1
            end

            states = Symbol[]
            state_values = [] 

            for nr in stable_roots
                H = 1.0 / (1.0 + (bi * nr)^c_p)
                dr = H / (αn * nr + bn)
                pheno = (nr * N0 > dr * D0) ? :R : :S
                push!(states, pheno)
                
                dim_N = nr * N0
                dim_D = dr * D0
                dim_I = (c_kT * dim_N * c_Dext) / c_gammaI
                push!(state_values, [dim_N, dim_D, dim_I])
            end

            if length(states) == 1
                label = states[1] == :S ? "Monostable {S}" : "Monostable {R}"
            else
                phenos_sorted = sort([string(s) for s in states])
                label = "Bistable {$(join(phenos_sorted, ","))}"
            end

            if label == "Bistable {S,S}"
                params_tuple = (c_betaN, c_betaD, c_gamma, c_gammaI, c_kC, c_kT, c_I0, c_Dext, c_Next, c_p)
                push!(ss_export_data, (p=params_tuple, states=state_values))
            end

            results[label] = get(results, label, 0) + 1
            
            if !haskey(env_counts, label)
                env_counts[label] = Dict(:Next_high => 0, :Dext_high => 0)
            end
            if c_Next > c_Dext
                env_counts[label][:Next_high] += 1
            else
                env_counts[label][:Dext_high] += 1
            end
        end
        next!(p_bar)
    end

    # --- 3. Output Results ---
    println("\n" * "="^85)
    println("GLOBAL ENVIRONMENT DISTRIBUTION (Across all $valid_count valid samples)")
    @printf("Next > Dext: %-10d (%.2f%%)\n", total_next_high, (total_next_high/valid_count)*100)
    @printf("Next < Dext: %-10d (%.2f%%)\n", total_dext_high, (total_dext_high/valid_count)*100)
    println("-"^85)
    @printf("%-30s | %-10s | %-18s | %-18s\n", "PHENOTYPE", "TOTAL %", "Next > Dext %", "Next < Dext %")
    println("-"^85)

    # Sort results by value in descending order
    sorted_keys = sort(collect(keys(results)), by=k->results[k], rev=true)

    for key in sorted_keys
        total_p = (results[key] / valid_count) * 100
        n_high_p = (env_counts[key][:Next_high] / results[key]) * 100
        d_high_p = (env_counts[key][:Dext_high] / results[key]) * 100
        @printf("%-30s | %-10.2f%% | %-18.2f%% | %-18.2f%%\n", key, total_p, n_high_p, d_high_p)
    end
    println("="^85)

    # Export details...
    open("SS_Phenotype_Parameters.txt", "w") do io
        println(io, "BISTABLE {S,S} PARAMETER AND STATE REPORT")
        println(io, "="^60)
        if isempty(ss_export_data)
            println(io, "No {S,S} cases found.")
        else
            for (idx, entry) in enumerate(ss_export_data)
                p = entry.p
                println(io, "CASE #$idx")
                @printf(io, "   Params: bN=%.2f, bD=%.2f, g=%.3f, gI=%.3f, kC=%.2e, kT=%.2e, I0=%.2f, Dext=%.2f, Next=%.2f, p=%.2f\n", 
                        p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9], p[10])
                for (s_idx, vals) in enumerate(entry.states)
                    @printf(io, "    State %d: N=%.2f, D=%.2f, I=%.2f\n", s_idx, vals[1], vals[2], vals[3])
                end
                println(io, "-"^40)
            end
        end
    end
    println("\n[Success] {S,S} details saved to 'SS_Phenotype_Parameters.txt'")
end

run_robust_parameter_sweep(Int(1e6))