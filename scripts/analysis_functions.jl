using Clustering, Statistics

function get_perfect_agent(environment::Environment)

    SD = environment[1][2]
    sigma = 2.0*log(SD) # convert SD to sigma

    agent = init_mixture_model(0.01, -6.0, sigma, 0.7)

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

function ambiguity_analysis(sigma::Float64 = 0.0, lambda::Float64 = 0.7, 
                           sd::Float64 = 1.0, n_timesteps::Int = 100; 
                           perfect_agent::Bool = false)
    DBMs = 0.5:0.5:6 # Distance between means
    n_iter = 10000

    mean_performances = Vector{Float64}(undef, length(DBMs))

    for (dbm_idx, DBM) in enumerate(DBMs)
        environment = get_simple_environment(DBM, sd)
        performances = Vector{Float64}(undef, n_iter)
        
        Threads.@threads for i in 1:n_iter
            if perfect_agent
                model = get_perfect_agent(environment)
            else
                model = simple_simulation(environment, sigma, lambda, n_timesteps)
            end
            performances[i] = get_behavioural_performance(model, environment)
        end
        
        mean_performances[dbm_idx] = mean(performances)
    end
    
    return mean_performances
end

function sigma_lambda_analysis(sigmas, lambda::Float64, 
                               distance_between_means::Float64, 
                               n_timesteps::Int = 100)
    n_iter = 10000
    environment = get_simple_environment(distance_between_means)
    
    mean_performances = Vector{Float64}(undef, length(sigmas))
    
    for (sigma_idx, sigma) in enumerate(sigmas)
        performances = Vector{Float64}(undef, n_iter)
        
        Threads.@threads for i in 1:n_iter
            model = simple_simulation(environment, sigma, lambda, n_timesteps)
            performances[i] = get_behavioural_performance(model, environment)
        end
        
        mean_performances[sigma_idx] = mean(performances)
    end
    
    return mean_performances
end


function get_trajectory(distance_between_means::Float64, change_points)
    environment = get_simple_environment(distance_between_means)
    model = simple_simulation(environment, 0.0, 0.7, 1)
    
    trajectory = Vector{Float64}(undef, 2000)
    trajectory[1] = get_behavioural_performance(model, environment)
    
    for t in 2:2000
        for cp in change_points
            if cp.timestep == t
                model.sigma = cp.sigma
                model.lambda = cp.lambda
            end
        end
        model = simulate_from_model(model, environment, 1)
        trajectory[t] = get_behavioural_performance(model, environment)
    end
    
    return trajectory
end

function trajectory_analysis(distance_between_means::Float64, timespan::Int, 
                            sigma::Float64 = 0.0, lambda::Float64 = 0.7)
    change_points = [
        (timestep = 200, sigma = sigma, lambda = lambda),
        (timestep = 200 + timespan, sigma = 0.0, lambda = 0.7)
    ]
    
    n_trajectories = 10000
    trajectories = Vector{Vector{Float64}}(undef, n_trajectories)
    
    Threads.@threads for i in 1:n_trajectories
        trajectories[i] = get_trajectory(distance_between_means, change_points)
    end
    
    return (change_points = change_points, trajectories = trajectories)
end

###############################################
## Helper function for outer parallelization ##
###############################################
function parallel_ambiguity_sigma(sigmas, sd::Float64, n_timesteps::Int)
    results = Vector{Vector{Float64}}(undef, length(sigmas))
    
    Threads.@threads for i in eachindex(sigmas)
        results[i] = ambiguity_analysis(sigmas[i], 0.7, sd, n_timesteps)
    end
    
    return results
end

function parallel_ambiguity_lambda(lambdas, sd::Float64, n_timesteps::Int)
    results = Vector{Vector{Float64}}(undef, length(lambdas))
    
    Threads.@threads for i in eachindex(lambdas)
        results[i] = ambiguity_analysis(0.0, lambdas[i], sd, n_timesteps)
    end
    
    return results
end

function parallel_sigma_lambda_analysis(sigmas, lambdas, dbm::Float64, n_timesteps::Int)
    results = Vector{Vector{Float64}}(undef, length(lambdas))
    
    Threads.@threads for i in eachindex(lambdas)
        results[i] = sigma_lambda_analysis(sigmas, lambdas[i], dbm, n_timesteps)
    end
    
    return results
end

function parallel_trajectory_sigma(sigmas, dbm::Float64, timespan::Int)
    results = Vector{Any}(undef, length(sigmas))
    
    Threads.@threads for i in eachindex(sigmas)
        results[i] = trajectory_analysis(dbm, timespan, sigmas[i], 0.7)
    end
    
    return results
end

function parallel_trajectory_lambda(lambdas, dbm::Float64, timespan::Int)
    results = Vector{Any}(undef, length(lambdas))
    
    Threads.@threads for i in eachindex(lambdas)
        results[i] = trajectory_analysis(dbm, timespan, 0.0, lambdas[i])
    end
    
    return results
end