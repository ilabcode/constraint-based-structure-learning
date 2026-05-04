const Environment = Vector{Tuple{Float64,Float64}}

function simple_simulation(environment::Environment, sigma::Float64 = 0.0, 
                          lambda::Float64 = 0.7, n_timesteps::Int = 100, g_map::Function = (x) -> 2.0 / pi * atan(x))
    model = init_mixture_model(.01, -6.0, sigma, lambda, g_map)
    simulate_from_model(model, environment, n_timesteps)
end


function get_simple_environment(distance_between_means::Float64, sd::Float64 = 1.0)
    DBM = distance_between_means
    Environment([(-DBM*sd, sd), (0,sd), (DBM*sd, sd)])
end

function simulate_from_model(model::MixModel, environment::Environment, n_timesteps::Int)
    
    environment_IDs = eachindex(environment)
    @inbounds for t in 1:n_timesteps
        i = rand(environment_IDs)
        observation = rand(Normal(environment[i]...))
        update_model!(model, observation)
    end

    return model
end

function simulate_for_animation(model_parameters, model_names, environment::Environment, n_timesteps)

    models = [init_mixture_model(para...) for para in model_parameters]

    update_sequence = Vector{Vector{MixModel}}(undef, n_timesteps)
    observations = Vector{Float64}(undef, n_timesteps)
    span = falses(n_timesteps)

    environment_IDs = eachindex(environment)

    for t in 1:n_timesteps
        
        i = rand(environment_IDs)
        observations[t] = rand(Normal(environment[i]...))

        for m in models
            update_model!(m, observations[t] )
        end
        
        update_sequence[t] = deepcopy(models)

    end

    return (update_sequence = update_sequence, 
            model_names = model_names, 
            observations = observations,
            span = span, 
            environment = environment)
end

function simulate_for_animation_w_param_change(span_parameters, model_names, environment::Environment, n_timesteps, span_length)

    models = [init_mixture_model(0.01,-6.0,0.0,0.7) for para in span_parameters]

    update_sequence = Vector{Vector{MixModel}}(undef, n_timesteps)
    observations = Vector{Float64}(undef, n_timesteps)
    span = falses(n_timesteps)
    span[200:(200+span_length)] .= true

    environment_IDs = eachindex(environment)

    for t in 1:n_timesteps
        
        i = rand(environment_IDs)
        observations[t] = rand(Normal(environment[i]...))

        for i in eachindex(models)
            if span[t]
                models[i].sigma = span_parameters[i].sigma
                models[i].lambda = span_parameters[i].lambda
            else
                models[i].sigma = 0.0
                models[i].lambda = 0.7
            end

            update_model!(models[i], observations[t])
        end
        
        update_sequence[t] = deepcopy(models)

    end

    return (update_sequence = update_sequence, 
            model_names = model_names, 
            observations = observations,
            span = span, 
            environment = environment)
end

function get_mixture_state(animation_data, name, timestep)
    model_index = findfirst(==(name), animation_data.model_names)
    animation_data.update_sequence[timestep][model_index].mixture
end


