# =====================================================================================
# GLOBAL ENVIRONMENT DISTRIBUTION (Across all 1000000 valid samples)
# Next > Dext: 110928     (11.09%)
# Next < Dext: 889072     (88.91%)
# -------------------------------------------------------------------------------------
# PHENOTYPE                      | TOTAL %    | Next > Dext %      | Next < Dext %
# -------------------------------------------------------------------------------------
# Bistable {R,S}                 | 45.47     % | 9.39              % | 90.61             %      
# Monostable {R}                 | 28.13     % | 7.48              % | 92.52             %      
# Monostable {S}                 | 25.48     % | 17.66             % | 82.34             %      
# Bistable {R,R}                 | 0.93      % | 23.79             % | 76.21             %      
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

    # All parameters now use linear uniform ranges
    uni_ranges = Dict(
        :betaN  => (250, 1000),
        :betaD  => (500,  2000),
        :kC     => (2.5e-4, 1e-3),
        :kT     => (2.5e-5, 1e-4),
        :Dext   => (500.0, 2000.0),
        :Next   => (250.0, 1000.0),
        :gamma  => (0.05, 0.2),
        :gammaI => (0.25, 1.0),
        :I0     => (100.0, 300.0),
        :p      => (1.0, 6.0)
    )

    p_bar = Progress(samples; dt=0.5, barglyphs=BarGlyphs("[=> ]"), barlen=40, color=:cyan)

    for i in 1:samples
        # Linear Uniform Sampling
        c_betaN  = rand() * (uni_ranges[:betaN][2] - uni_ranges[:betaN][1]) + uni_ranges[:betaN][1]
        c_betaD  = rand() * (uni_ranges[:betaD][2] - uni_ranges[:betaD][1]) + uni_ranges[:betaD][1]
        c_kC     = rand() * (uni_ranges[:kC][2] - uni_ranges[:kC][1]) + uni_ranges[:kC][1]
        c_kT     = rand() * (uni_ranges[:kT][2] - uni_ranges[:kT][1]) + uni_ranges[:kT][1]
        c_Dext   = rand() * (uni_ranges[:Dext][2] - uni_ranges[:Dext][1]) + uni_ranges[:Dext][1]
        c_Next   = rand() * (uni_ranges[:Next][2] - uni_ranges[:Next][1]) + uni_ranges[:Next][1]
        
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

    sorted_keys = sort(collect(keys(results)), by=k->results[k], rev=true)

    for key in sorted_keys
        total_p = (results[key] / valid_count) * 100
        n_high_p = (env_counts[key][:Next_high] / results[key]) * 100
        d_high_p = (env_counts[key][:Dext_high] / results[key]) * 100
        @printf("%-30s | %-10.2f%% | %-18.2f%% | %-18.2f%%\n", key, total_p, n_high_p, d_high_p)
    end
    println("="^85)

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