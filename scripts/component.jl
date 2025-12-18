######################
## COMPONENT STRUCT ##
######################
struct ParentNode
    μ::Float64
    μ̂::Float64
    π::Float64
    π̂::Float64
    ω::Float64
end

mutable struct Component
    value_parent::ParentNode
    noise_parent::ParentNode
    zeta::Float64
    weight::Float64
    theta::Tuple{Float64, Float64}  # (mean, sd)
    new::Bool
end

function new_component(mean::Float64, sigma::Float64, omega::Float64, alpha::Float64)
    value_parent = ParentNode(mean, 0.0, 1.0, 0.0, omega)
    noise_parent = ParentNode(0.0, 0.0, 1.0, 0.0, omega)
    
    weight = alpha
    component = Component(value_parent, noise_parent, sigma, weight, (0.0,0.0), true)
    component.theta = get_component_theta(component)
    
    return component
end

#######################
## UTILITY FUNCTIONS ##
#######################

function get_component_theta(component::Component)

    mean = component.value_parent.μ
    precision = 1.0 / exp(component.zeta + component.noise_parent.μ)
    sd = 1.0 / sqrt(precision)
    
    return (mean, sd)
end

######################
## UPDATE FUNCTIONS ##
######################
function update_prediction(parent::ParentNode) 
    μ̂ = parent.μ
    π̂ = 1.0 / (1.0 / parent.π + exp(parent.ω))
    return ParentNode(parent.μ, μ̂, parent.π, π̂, parent.ω)
end

function update_component!(component::Component, observation::Float64, activity_level::Float64, lambda::Float64)

    ζ = component.zeta
    λ = lambda
    a = activity_level
    
    if component.new
        component.weight = a
        component.new = false
    else
        ### Component level ###
        vp = update_prediction(component.value_parent)
        np = update_prediction(component.noise_parent)

        π̂ᴵ = 1.0 / exp(ζ + np.μ)

        δ = observation - vp.μ̂ # Value prediction error
        ϵ = π̂ᴵ/vp.π + π̂ᴵ * δ^2 - 1.0 # Noise prediction error

        # Value parent updates
        vp_π = vp.π̂ + π̂ᴵ * a
        vp_μ = vp.μ̂ + δ * π̂ᴵ / vp_π * a
        
        # Noise parent updates  
        np_π = np.π̂ + 0.5 * (1.0 + ϵ) * a
        np_μ = np.μ̂ + 0.5 * (1.0 / np_π) * ϵ * a

        # Store updated values
        component.value_parent = ParentNode(vp_μ, vp.μ̂, vp_π, vp.π̂, vp.ω)
        component.noise_parent = ParentNode(np_μ, np.μ̂, np_π, np.π̂, np.ω)

        ### System level ###
        θ = get_component_theta(component)
        component.theta = θ

        τ = 1.0 / (θ[2]^2) # converting from SD to precision
        μ = θ[1]
        g = 2.0 / pi * atan(exp(λ) / ((observation-μ)^2 * τ))
        component.weight += g * (a - component.weight)
    end

    return nothing
end



