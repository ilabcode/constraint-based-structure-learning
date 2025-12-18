using Distributions, StatsFuns

#############
## STRUCTS ##
#############
mutable struct GaussianMixture
    components::Vector{Component}
    IDs::Vector{Int}
    next_ID::Int
end

function init_GaussianMixture()
    GaussianMixture(Component[], Int[], 1)
end

mutable struct MixModel
    mixture::GaussianMixture
    alpha::Float64
    omega::Float64
    sigma::Float64
    lambda::Float64
end

function init_mixture_model(alpha::Float64 = 0.01, omega::Float64 = -6.0, 
                           sigma::Float64 = 0.0, lambda::Float64 = 0.7)
    MixModel(init_GaussianMixture(), alpha, omega, sigma, lambda)
end

#######################
## UTILITY FUNCTIONS ##
#######################
function add_component_to_mixture!(model::MixModel, mean::Float64)
    component = new_component(mean, model.sigma, model.omega, model.alpha)
    push!(model.mixture.components, component)
    push!(model.mixture.IDs, model.mixture.next_ID)
    model.mixture.next_ID += 1
    return nothing
end

function get_activity_levels(model::MixModel, observation::Float64)
    activity_levels = Vector{Float64}(undef, length(model.mixture.components))
    
    o = observation
    components = model.mixture.components
    normalisation = 0.0

    @inbounds for i in eachindex(components)
        c = components[i]
        if c.new
            val = model.alpha * pdf(Normal(c.theta[1], c.theta[2]), o) * c.weight
        else
            τ = 1.0 / (c.theta[2]^2)
            g = 2.0 / π * atan(exp(model.lambda) / ((o - c.theta[1])^2 * τ))
            f = pdf(Normal(c.theta[1], c.theta[2]), o) 
            val = g * f * c.weight
        end
        activity_levels[i] = val
        normalisation += val
    end
    
    # Normalise
    inv_norm = 1.0 / normalisation
    @inbounds for i in eachindex(components)
        activity_levels[i] *= inv_norm
    end

    return activity_levels
end


######################
## UPDATE FUNCTIONS ##
######################
function update_model!(model, observation)

    add_component_to_mixture!(model, observation)
    
    activity_levels = get_activity_levels(model, observation)

    components = model.mixture.components

    # Update all components
    for i in eachindex(components)
        update_component!(
            components[i], 
            observation, 
            activity_levels[i], 
            model.lambda
            )
    end

    # Filter out dead components. Keep surviving components
    suvivors = [c.weight >= model.alpha for c in components]
    if !all(suvivors)
        filter = findall(suvivors)
        model.mixture.components = components[filter]
        model.mixture.IDs = model.mixture.IDs[filter]
    end

    return nothing

end


