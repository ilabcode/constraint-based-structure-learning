include("../component.jl")
include("../mixture.jl")
include("../simulation_functions.jl")
include("../analysis_functions.jl")
include("../export_functions.jl")
include("../plot_functions.jl")

g_map = (x) -> tanh(x / 2.0)
data_dir = joinpath(@__DIR__, "../../data")
suffix = "logistic"

##########################
### AMBIGUITY ANALYSIS ###
##########################
ambiguity_sigmas = -2.2:.2:2
ambiguity_lambdas = -2.5:.3:4
n_timesteps_short = 100
n_timesteps_long = 10000
n_iter = 1000
baseline_sigma = 0.0
baseline_lambda = 0.7

perfect_agent_ambiguity_results = ambiguity_analysis(perfect_agent = true, g_map = g_map, n_iter = n_iter)
CSV.write(joinpath(data_dir, "perfect_agent_results_$(suffix).csv"), ambiguity_to_dataframe([perfect_agent_ambiguity_results]))
println("perfect_agent_results_$(suffix).csv saved.")

convergence_baseline_results = ambiguity_analysis(0.0, 0.7, 1.0, 100000, g_map = g_map, n_iter = n_iter)
CSV.write(joinpath(data_dir, "convergence_baseline_results_$(suffix).csv"), ambiguity_to_dataframe([convergence_baseline_results]))
println("convergence_baseline_results_$(suffix).csv saved.")

# Short term
ambiguity_short_sd1 = parallel_ambiguity_analysis(ambiguity_sigmas, baseline_sigma, ambiguity_lambdas, baseline_lambda, 1.0, n_timesteps_short, g_map, n_iter=n_iter)

# Long term
ambiguity_long_sd1 = parallel_ambiguity_analysis(ambiguity_sigmas, baseline_sigma, ambiguity_lambdas, baseline_lambda, 1.0, n_timesteps_long, g_map, n_iter=n_iter)

# Save ambiguity results
CSV.write(joinpath(data_dir, "ambiguity_results_$(suffix).csv"), ambiguity_to_dataframe(vcat(
    ambiguity_short_sd1, ambiguity_long_sd1
)))
println("ambiguity_results_$(suffix).csv saved.")


#########################################
### SIGMA LAMBDA INTERACTION ANALYSIS ###
#########################################
interaction_sigmas = -4.2:.2:3.6
interaction_lambdas = -3.1:.2:4.7

n_iter_sigma_lambda = n_iter

# Unambiguous
sigma_lambda_short_DBM6 = parallel_sigma_lambda_analysis(interaction_sigmas, interaction_lambdas, 6.0, n_timesteps_short, g_map, n_iter=n_iter_sigma_lambda)
sigma_lambda_long_DBM6  = parallel_sigma_lambda_analysis(interaction_sigmas, interaction_lambdas, 6.0, n_timesteps_long, g_map, n_iter=n_iter_sigma_lambda)

# Semi-ambiguous
sigma_lambda_short_DBM45 = parallel_sigma_lambda_analysis(interaction_sigmas, interaction_lambdas, 4.5, n_timesteps_short, g_map, n_iter=n_iter_sigma_lambda)
sigma_lambda_long_DBM45  = parallel_sigma_lambda_analysis(interaction_sigmas, interaction_lambdas, 4.5, n_timesteps_long, g_map, n_iter=n_iter_sigma_lambda)

# Ambiguous
sigma_lambda_short_DBM3 = parallel_sigma_lambda_analysis(interaction_sigmas, interaction_lambdas, 3.0, n_timesteps_short, g_map, n_iter=n_iter_sigma_lambda)
sigma_lambda_long_DBM3  = parallel_sigma_lambda_analysis(interaction_sigmas, interaction_lambdas, 3.0, n_timesteps_long, g_map, n_iter=n_iter_sigma_lambda)

# Save interaction results
CSV.write(joinpath(data_dir, "sigma_lambda_results_$(suffix).csv"), sigma_lambda_to_dataframe(vcat(
    sigma_lambda_short_DBM6, sigma_lambda_long_DBM6,
    sigma_lambda_short_DBM45, sigma_lambda_long_DBM45,
    sigma_lambda_short_DBM3, sigma_lambda_long_DBM3
)))
println("sigma_lambda_results_$(suffix).csv saved.")


###########################
### TRAJECTORY ANALYSIS ###
###########################
trajectory_sigmas = 0.0:-0.8:-4
trajectory_lambdas = 0.7:-0.6:-2.3

n_iter_trajectory = n_iter

# Unambiguous
trajectory_DBM6_span25  = parallel_trajectory_analysis(trajectory_sigmas, baseline_sigma, trajectory_lambdas, baseline_lambda, 6.0, 25, n_iter=n_iter_trajectory, g_map=g_map)
trajectory_DBM6_span100 = parallel_trajectory_analysis(trajectory_sigmas, baseline_sigma, trajectory_lambdas, baseline_lambda, 6.0, 100, n_iter=n_iter_trajectory, g_map=g_map)
trajectory_DBM6_span400 = parallel_trajectory_analysis(trajectory_sigmas, baseline_sigma, trajectory_lambdas, baseline_lambda, 6.0, 400, n_iter=n_iter_trajectory, g_map=g_map)

# Ambiguous
trajectory_DBM3_span25  = parallel_trajectory_analysis(trajectory_sigmas, baseline_sigma, trajectory_lambdas, baseline_lambda, 3.0, 25, n_iter=n_iter_trajectory, g_map=g_map)
trajectory_DBM3_span100 = parallel_trajectory_analysis(trajectory_sigmas, baseline_sigma, trajectory_lambdas, baseline_lambda, 3.0, 100, n_iter=n_iter_trajectory, g_map=g_map)
trajectory_DBM3_span400 = parallel_trajectory_analysis(trajectory_sigmas, baseline_sigma, trajectory_lambdas, baseline_lambda, 3.0, 400, n_iter=n_iter_trajectory, g_map=g_map)

# Save trajectory results
CSV.write(joinpath(data_dir, "trajectory_results_$(suffix).csv"), trajectory_to_dataframe(vcat(
    trajectory_DBM6_span25, trajectory_DBM6_span100, trajectory_DBM6_span400,
    trajectory_DBM3_span25, trajectory_DBM3_span100, trajectory_DBM3_span400
)))
println("trajectory_results_$(suffix).csv saved.")
