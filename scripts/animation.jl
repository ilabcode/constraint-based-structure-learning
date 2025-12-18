include("simulation_functions.jl")
include("plot_functions.jl")

# Setup
unambiguous = environment = get_simple_environment(6.0)
ambiguous = environment = get_simple_environment(3.0)

standard_alpha = 0.01
standard_omega = -6.0
standard_sigma = 0.0
standard_lambda = 0.7

# Sigma
sigma_names = ["σ = 2.0","σ = 0.5", "σ = 0.0", "σ = -0.5", "σ = -1.0", "σ = -2.0"]

sigma_params = [
    (standard_alpha, standard_omega, 2.0, standard_lambda),
    (standard_alpha, standard_omega, 0.5, standard_lambda),
    (standard_alpha, standard_omega, 0.0, standard_lambda),
    (standard_alpha, standard_omega, -0.5, standard_lambda),
    (standard_alpha, standard_omega, -1.0, standard_lambda),
    (standard_alpha, standard_omega, -2.0, standard_lambda)
]

animate_simulation(
    simulate_for_animation(
        sigma_params, sigma_names, 
        unambiguous, 200
        ); 
    show_data=true
    )

animate_simulation(
    simulate_for_animation(
        sigma_params, sigma_names, 
        ambiguous, 200
        ); 
    show_data=true
    )



# Lambda
lambda_names = ["λ = 4.0","λ = 2.0", "λ = 0.7", "λ = 0.0", "λ = -1.0", "λ = -2.0"]

lambda_params = [
    (standard_alpha, standard_omega, standard_sigma, 4.0),
    (standard_alpha, standard_omega, standard_sigma, 2.0),
    (standard_alpha, standard_omega, standard_sigma, 0.7),
    (standard_alpha, standard_omega, standard_sigma, 0.0),
    (standard_alpha, standard_omega, standard_sigma, -1.0),
    (standard_alpha, standard_omega, standard_sigma, -2.0)
]

animate_simulation(
    simulate_for_animation(
        lambda_params, lambda_names, 
        unambiguous, 200
        ); 
    show_data=true
    )

animate_simulation(
    simulate_for_animation(
        lambda_params, lambda_names, 
        ambiguous, 200
        ); 
    show_data=true
    )


# Long term
long_names = ["standard","σ = -1.0", "σ = -2.0", "λ = 0.0", "λ = -1.0", "λ = -2.0"]

long_params = [
    (standard_alpha, standard_omega, standard_sigma, standard_lambda),
    (standard_alpha, standard_omega, -1.0, standard_lambda),
    (standard_alpha, standard_omega, -2.0, standard_lambda),
    (standard_alpha, standard_omega, standard_sigma, 0.0),
    (standard_alpha, standard_omega, standard_sigma, -1.0),
    (standard_alpha, standard_omega, standard_sigma, -2.0)
]

animate_simulation(
    simulate_for_animation(
        long_params, long_names, 
        unambiguous, 2000
        ); 
    show_data=false,
    speed = 10
    )

animate_simulation(
    simulate_for_animation(
        long_params, long_names, 
        ambiguous, 2000
        ); 
    show_data=false,
    speed = 10
    )



# Trejectory
trejactory_names = ["standard","σ = -2.4", "σ = -4.0", "λ = 0.1", "λ = -0.5", "λ = -1.7"]

trejactory_params = [
    (sigma = standard_sigma, lambda = standard_lambda),
    (sigma = -2.4, lambda = standard_lambda),
    (sigma = -4.0, lambda = standard_lambda),
    (sigma = standard_sigma, lambda = 0.1),
    (sigma = standard_sigma, lambda = -0.5),
    (sigma = standard_sigma, lambda = -1.7)
]


animate_simulation(
    simulate_for_animation_w_param_change(
        trejactory_params, trejactory_names, 
        unambiguous, 1000, 25
        );
    show_data=false,
    speed = 5
    )

animate_simulation(
    simulate_for_animation_w_param_change(
        trejactory_params, trejactory_names, 
        unambiguous, 1000, 100
        );
    show_data=false,
    speed = 5
    )

animate_simulation(
    simulate_for_animation_w_param_change(
        trejactory_params, trejactory_names, 
        unambiguous, 1000, 400
        );
    show_data=false,
    speed = 5
    )

animate_simulation(
    simulate_for_animation_w_param_change(
        trejactory_params, trejactory_names, 
        ambiguous, 1000, 25
        );
    show_data=false,
    speed = 5
    )

animate_simulation(
    simulate_for_animation_w_param_change(
        trejactory_params, trejactory_names, 
        ambiguous, 1000, 100
        );
    show_data=false,
    speed = 5
    )

animate_simulation(
    simulate_for_animation_w_param_change(
        trejactory_params, trejactory_names, 
        ambiguous, 1000, 400
        );
    show_data=false,
    speed = 5
    )

