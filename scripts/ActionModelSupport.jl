
using ActionModels, Distributions

# ========================
# SUBMODEL STRUCT
# ========================

Base.@kwdef struct MixModelSubmodel <: ActionModels.AbstractSubmodel
    sigma::Float64  = 0.0
    lambda::Float64 = 0.7
end

function ActionModels.get_parameter_types(::MixModelSubmodel)
    return (sigma = Float64, lambda = Float64)
end

function ActionModels.get_state_types(::MixModelSubmodel)
    return (n_components = Int64,)
end


# ========================
# ATTRIBUTES STRUCT
# ========================

mutable struct MixModelAttributes{T<:Real} <: ActionModels.AbstractSubmodelAttributes
    sigma::T
    lambda::T
    mixture::GaussianMixture
    n_components::Int
end

function ActionModels.initialize_attributes(
    submodel::MixModelSubmodel,
    ::Type{TF} = Float64,
    ::Type{TI} = Int64,
) where {TF, TI}
    MixModelAttributes{TF}(
        TF(submodel.sigma),
        TF(submodel.lambda),
        init_GaussianMixture(),
        0,
    )
end


# ========================
# REQUIRED API METHODS
# ========================

function ActionModels.reset!(attributes::MixModelAttributes)
    attributes.mixture = init_GaussianMixture()
    attributes.n_components = 0
    return nothing
end

function ActionModels.get_parameters(attributes::MixModelAttributes)
    return (; sigma = attributes.sigma, lambda = attributes.lambda)
end

function ActionModels.get_states(attributes::MixModelAttributes)
    return (; n_components = attributes.n_components)
end

function ActionModels.get_parameters(attributes::MixModelAttributes, name::Symbol)
    name in (:sigma, :lambda) ? getfield(attributes, name) : AttributeError()
end

function ActionModels.get_states(attributes::MixModelAttributes, name::Symbol)
    name == :n_components ? attributes.n_components : AttributeError()
end

function ActionModels.set_parameters!(attributes::MixModelAttributes, name::Symbol, value::T) where {T<:Real}
    name in (:sigma, :lambda) ? (setfield!(attributes, name, value); true) : AttributeError()
end

function ActionModels.set_states!(attributes::MixModelAttributes, name::Symbol, value::T) where {T<:Real}
    name == :n_components ? (attributes.n_components = Int(round(Float64(value))); true) : AttributeError()
end


# ========================
# UPDATE
# ========================

function update!(attributes::MixModelAttributes, observation::Float64, sigma::Float64, lambda::Float64)
    g_map = x -> 2.0 / pi * atan(x)
    model = MixModel(
        attributes.mixture,
        0.01,
        -6.0,
        sigma,
        lambda,
        g_map
    )
    activity_levels = update_model!(model, observation)
    attributes.mixture    = model.mixture
    attributes.n_components = length(model.mixture.components)
    return (activity_levels, model.mixture.IDs)
end


# ========================
# ACTION MODEL FUNCTION
# ========================

function constraint_based_structure_learning(attributes::ModelAttributes, observation::Float64)
    mix = attributes.submodel
    sigma_val  = Float64(attributes.parameters.sigma.value)
    lambda_val = Float64(attributes.parameters.lambda.value)

    activity_levels, IDs = update!(mix, observation, sigma_val, lambda_val)

    # n_timesteps = 100

    # categorical_dist = fill(eps(), n_timesteps)
    # for i in eachindex(IDs)
    #     if isfinite(activity_levels[i])
    #         categorical_dist[IDs[i]] += activity_levels[i]
    #     end
    # end

    # categorical_dist ./= sum(categorical_dist)

    action_distribution = Categorical(activity_levels)

    return action_distribution
end
