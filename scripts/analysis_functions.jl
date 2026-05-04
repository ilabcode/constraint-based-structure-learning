using Clustering, Statistics, ExactOptimalTransport, LinearAlgebra, HiGHS, MathOptInterface, CSV
const MOI = MathOptInterface

function MW2(mix1::GaussianMixture, mix2::GaussianMixture)
    K1 = length(mix1.components)
    K2 = length(mix2.components)

    # Normalize weights to simplex
    w1 = [c.weight for c in mix1.components]
    w2 = [c.weight for c in mix2.components]
    π1 = w1 ./ sum(w1)
    π2 = w2 ./ sum(w2)

    # W₂² cost matrix — 1D closed form: (Δm)² + (Δσ)²
    C = Matrix{Float64}(undef, K1, K2)
    for k in 1:K1
        m1, σ1 = mix1.components[k].theta
        for l in 1:K2
            m2, σ2 = mix2.components[l].theta
            C[k, l] = (m1 - m2)^2 + (σ1 - σ2)^2
        end
    end

    # Earth Mover's Distance (W₂²) via linear programming
    solver = MOI.instantiate(MOI.OptimizerWithAttributes(HiGHS.Optimizer, MOI.Silent() => true))
    γ = emd(π1, π2, C, solver)
    return sqrt(dot(γ, C))
end

function get_perfect_agent(environment::Environment, g_map::Function = (x) -> 2.0 / pi * atan(x))

    SD = environment[1][2]
    sigma = 2.0*log(SD) # convert SD to sigma

    agent = init_mixture_model(0.01, -6.0, sigma, 0.7, g_map)

    for theta in environment
        add_component_to_mixture!(agent, theta[1])
        agent.mixture.components[end].weight = 1.0
        agent.mixture.components[end].new = false
    end

    return agent

end

function get_behavioural_performance(model::MixModel, environment::Environment)
    # Pre-allocate arrays
    obs_per_env_comp = 20
    n_total = length(environment) * obs_per_env_comp
    actions = Vector{Int}(undef, n_total)
    ground_truths = Vector{Int}(undef, n_total)

    idx = 1
    @inbounds for env_i in eachindex(environment)
        mean, sd = environment[env_i]
        normal_dist = Normal(mean, sd)
        
        for i in 1:obs_per_env_comp
            observation = rand(normal_dist)
            activity_levels = get_activity_levels(model, observation)

            if any(isnan, activity_levels) 
                # Using uniform categorical distribution
                # when activity levels fail due to calculation without running the full algorithm
                # Happens when an observation falls very far from all components (and here no new components are added)
                # Usually this is due to extremly low sigma
                action_distribution = eachindex(activity_levels)
            else
                # The distribution used in all normal cases
                action_distribution = Categorical(activity_levels)
            end

            action = rand(action_distribution)
            
            actions[idx] = action
            ground_truths[idx] = env_i
            idx += 1
        end
    end
    
    return randindex(actions, ground_truths)[1]
end

function summarise(v::Vector{Float64})
    sorted = sort(v)
    n = length(sorted)
    (mean   = mean(v),
     min    = sorted[1],
     q01    = sorted[max(1, round(Int, 0.01n))],
     q05    = sorted[max(1, round(Int, 0.05n))],
     q25    = sorted[max(1, round(Int, 0.25n))],
     median = sorted[max(1, round(Int, 0.50n))],
     q75    = sorted[max(1, round(Int, 0.75n))],
     q95    = sorted[max(1, round(Int, 0.95n))],
     q99    = sorted[max(1, round(Int, 0.99n))],
     max    = sorted[n])
end

function ambiguity_analysis(sigma::Float64 = 0.0, lambda::Float64 = 0.7,
                           sd::Float64 = 1.0, n_timesteps::Int = 100;
                           perfect_agent::Bool = false, g_map::Function = (x) -> 2.0 / pi * atan(x),
                           n_iter::Int = 10000)
    DBMs = 0.5:0.5:6 # Distance between means
    n_dbm = length(DBMs)

    ari_stats   = Vector{NamedTuple}(undef, n_dbm)
    mw2_stats   = Vector{NamedTuple}(undef, n_dbm)
    ncomp_stats = Vector{NamedTuple}(undef, n_dbm)

    for (dbm_idx, DBM) in enumerate(DBMs)
        environment = get_simple_environment(DBM, sd)
        perfect = get_perfect_agent(environment, g_map)

        ari_vals   = Vector{Float64}(undef, n_iter)
        mw2_vals   = Vector{Float64}(undef, n_iter)
        ncomp_vals = Vector{Float64}(undef, n_iter)

        Threads.@threads for i in 1:n_iter
            if perfect_agent
                model = perfect
            else
                model = simple_simulation(environment, sigma, lambda, n_timesteps, g_map)
            end
            ari_vals[i]   = get_behavioural_performance(model, environment)
            mw2_vals[i]   = MW2(model.mixture, perfect.mixture)
            ncomp_vals[i] = Float64(length(model.mixture.components))
        end

        ari_stats[dbm_idx]   = summarise(ari_vals)
        mw2_stats[dbm_idx]   = summarise(mw2_vals)
        ncomp_stats[dbm_idx] = summarise(ncomp_vals)
    end

    return (DBMs = collect(DBMs), ari = ari_stats, mw2 = mw2_stats, ncomp = ncomp_stats, sigma = sigma, lambda = lambda, sd = sd, n_timesteps = n_timesteps)
end

function sigma_lambda_analysis(sigmas, lambda::Float64,
                               distance_between_means::Float64,
                               n_timesteps::Int = 100;
                               g_map::Function = (x) -> 2.0 / pi * atan(x),
                               n_iter::Int = 10000)
    environment = get_simple_environment(distance_between_means)
    perfect = get_perfect_agent(environment, g_map)
    n_sigmas = length(sigmas)

    ari_stats   = Vector{NamedTuple}(undef, n_sigmas)
    mw2_stats   = Vector{NamedTuple}(undef, n_sigmas)
    ncomp_stats = Vector{NamedTuple}(undef, n_sigmas)

    for (sigma_idx, sigma) in enumerate(sigmas)
        ari_vals   = Vector{Float64}(undef, n_iter)
        mw2_vals   = Vector{Float64}(undef, n_iter)
        ncomp_vals = Vector{Float64}(undef, n_iter)

        Threads.@threads for i in 1:n_iter
            model = simple_simulation(environment, sigma, lambda, n_timesteps, g_map)
            ari_vals[i]   = get_behavioural_performance(model, environment)
            mw2_vals[i]   = MW2(model.mixture, perfect.mixture)
            ncomp_vals[i] = Float64(length(model.mixture.components))
        end

        ari_stats[sigma_idx]   = summarise(ari_vals)
        mw2_stats[sigma_idx]   = summarise(mw2_vals)
        ncomp_stats[sigma_idx] = summarise(ncomp_vals)
    end

    return (sigmas = collect(sigmas), ari = ari_stats, mw2 = mw2_stats, ncomp = ncomp_stats, lambda = lambda, distance_between_means = distance_between_means, n_timesteps = n_timesteps)
end


function get_trajectory(distance_between_means::Float64, change_points, perfect::MixModel; g_map::Function = (x) -> 2.0 / pi * atan(x))
    environment = get_simple_environment(distance_between_means)
    model = simple_simulation(environment, 0.0, 0.7, 1, g_map)

    T = 2000
    ari   = Vector{Float64}(undef, T)
    mw2   = Vector{Float64}(undef, T)
    ncomp = Vector{Float64}(undef, T)

    ari[1]   = get_behavioural_performance(model, environment)
    mw2[1]   = MW2(model.mixture, perfect.mixture)
    ncomp[1] = Float64(length(model.mixture.components))

    for t in 2:T
        for cp in change_points
            if cp.timestep == t
                model.sigma = cp.sigma
                model.lambda = cp.lambda
            end
        end
        model = simulate_from_model(model, environment, 1)
        ari[t]   = get_behavioural_performance(model, environment)
        mw2[t]   = MW2(model.mixture, perfect.mixture)
        ncomp[t] = Float64(length(model.mixture.components))
    end

    return (ari = ari, mw2 = mw2, ncomp = ncomp)
end

function trajectory_analysis(distance_between_means::Float64, timespan::Int,
                            sigma::Float64 = 0.0, lambda::Float64 = 0.7;
                            n_iter::Int = 10000, g_map::Function = (x) -> 2.0 / pi * atan(x))
    change_points = [
        (timestep = 200, sigma = sigma, lambda = lambda),
        (timestep = 200 + timespan, sigma = 0.0, lambda = 0.7)
    ]

    environment = get_simple_environment(distance_between_means)
    perfect = get_perfect_agent(environment, g_map)
    T = 2000

    ari_mat   = Matrix{Float64}(undef, T, n_iter)
    mw2_mat   = Matrix{Float64}(undef, T, n_iter)
    ncomp_mat = Matrix{Float64}(undef, T, n_iter)

    Threads.@threads for i in 1:n_iter
        traj = get_trajectory(distance_between_means, change_points, perfect, g_map = g_map)
        ari_mat[:, i]   = traj.ari
        mw2_mat[:, i]   = traj.mw2
        ncomp_mat[:, i] = traj.ncomp
    end

    ari_stats   = [summarise(Vector(ari_mat[t, :]))   for t in 1:T]
    mw2_stats   = [summarise(Vector(mw2_mat[t, :]))   for t in 1:T]
    ncomp_stats = [summarise(Vector(ncomp_mat[t, :]))  for t in 1:T]

    return (change_points = change_points, ari = ari_stats, mw2 = mw2_stats, ncomp = ncomp_stats,
            distance_between_means = distance_between_means, sigma = sigma, lambda = lambda)
end

function convergence_investigation(distance_between_means::Float64; n_iter::Int = 10000, g_map::Function = (x) -> 2.0 / pi * atan(x))
    environment = get_simple_environment(distance_between_means)

    all_performance  = Matrix{Float64}(undef, 1000, n_iter)
    all_n_components = Matrix{Float64}(undef, 1000, n_iter)
    all_sum_weights  = Matrix{Float64}(undef, 1000, n_iter)

    Threads.@threads for j in 1:n_iter
        model = simple_simulation(environment, 0.0, 0.7, 100, g_map)
        for i in 1:1000
            model = simulate_from_model(model, environment, 10)
            all_performance[i, j]  = get_behavioural_performance(model, environment)
            all_n_components[i, j] = length(model.mixture.components)
            all_sum_weights[i, j]  = sum(c.weight for c in model.mixture.components)
        end
    end

    return (performance  = vec(mean(all_performance,  dims=2)),
            n_components = vec(mean(all_n_components, dims=2)),
            sum_weights  = vec(mean(all_sum_weights,  dims=2)))
end


###############################################
## Helper function for outer parallelization ##
###############################################
function parallel_ambiguity_analysis(sigmas, baseline_sigma::Float64, lambdas, baseline_lambda::Float64,
                                      sd::Float64, n_timesteps::Int, g_map::Function = (x) -> 2.0 / pi * atan(x); n_iter::Int = 10000)
    # Sigma sweep (all sigmas at baseline lambda)
    sigma_results = Vector{Any}(undef, length(sigmas))
    Threads.@threads for i in eachindex(sigmas)
        sigma_results[i] = ambiguity_analysis(sigmas[i], baseline_lambda, sd, n_timesteps, g_map = g_map, n_iter = n_iter)
    end

    # Lambda sweep (all lambdas at baseline sigma, skip the intersection)
    lambdas_filtered = filter(l -> l != baseline_lambda, lambdas)
    lambda_results = Vector{Any}(undef, length(lambdas_filtered))
    Threads.@threads for i in eachindex(lambdas_filtered)
        lambda_results[i] = ambiguity_analysis(baseline_sigma, lambdas_filtered[i], sd, n_timesteps, g_map = g_map, n_iter = n_iter)
    end

    return vcat(sigma_results, lambda_results)
end

function parallel_sigma_lambda_analysis(sigmas, lambdas, dbm::Float64, n_timesteps::Int, g_map::Function = (x) -> 2.0 / pi * atan(x); n_iter::Int = 10000)
    results = Vector{Any}(undef, length(lambdas))

    Threads.@threads for i in eachindex(lambdas)
        results[i] = sigma_lambda_analysis(sigmas, lambdas[i], dbm, n_timesteps, g_map = g_map, n_iter = n_iter)
    end

    return results
end

function parallel_trajectory_analysis(sigmas, baseline_sigma::Float64, lambdas, baseline_lambda::Float64,
                                       dbm::Float64, timespan::Int; n_iter::Int = 10000, g_map::Function = (x) -> 2.0 / pi * atan(x))
    # Sigma sweep (all sigmas at baseline lambda)
    sigma_results = Vector{Any}(undef, length(sigmas))
    Threads.@threads for i in eachindex(sigmas)
        sigma_results[i] = trajectory_analysis(dbm, timespan, sigmas[i], baseline_lambda, n_iter = n_iter, g_map = g_map)
    end

    # Lambda sweep (all lambdas at baseline sigma, skip the intersection)
    lambdas_filtered = filter(l -> l != baseline_lambda, lambdas)
    lambda_results = Vector{Any}(undef, length(lambdas_filtered))
    Threads.@threads for i in eachindex(lambdas_filtered)
        lambda_results[i] = trajectory_analysis(dbm, timespan, baseline_sigma, lambdas_filtered[i], n_iter = n_iter, g_map = g_map)
    end

    return vcat(sigma_results, lambda_results)
end