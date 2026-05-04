include("component.jl")
include("mixture.jl")
include("simulation_functions.jl")
include("analysis_functions.jl")
include("export_functions.jl")
include("plot_functions.jl")

# G map Functions
arctan = (x) -> 2.0 / pi * atan(x)
logistic = (x) -> tanh(x/2.0)
exp_decay = (x) -> 1-exp(-x)
#Set function here:
g_map = arctan


test = convergence_investigation(6.0, n_iter = 10000, g_map = g_map)
cor(test[3], test[2])
plot(test[1])
plot(test[2])
plot(test[3])

##########################
### AMBIGUITY ANALYSIS ###
##########################
ambiguity_sigmas = -2.2:.2:2
ambiguity_lambdas = -2.5:.3:4 
n_timesteps_short = 100
n_timesteps_long = 10000
n_iter = 10000
baseline_sigma = 0.0
baseline_lambda = 0.7

perfect_agent_ambiguity_results = ambiguity_analysis(perfect_agent = true, g_map = g_map, n_iter = n_iter)
CSV.write("data/perfect_agent_results.csv", ambiguity_to_dataframe([perfect_agent_ambiguity_results]))
println("perfect_agent_results.csv saved.")

convergence_baseline_results = ambiguity_analysis(0.0, 0.7, 1.0, 100000, g_map = g_map, n_iter = n_iter)
CSV.write("data/convergence_baseline_results.csv", ambiguity_to_dataframe([convergence_baseline_results]))
println("convergence_baseline_results.csv saved.")

# Short term
ambiguity_short_sd1  = parallel_ambiguity_analysis(ambiguity_sigmas, baseline_sigma, ambiguity_lambdas, baseline_lambda, 1.0, n_timesteps_short, g_map, n_iter=n_iter)
ambiguity_short_sd2  = parallel_ambiguity_analysis(ambiguity_sigmas, baseline_sigma, ambiguity_lambdas, baseline_lambda, 2.0, n_timesteps_short, g_map, n_iter=n_iter)
ambiguity_short_sd05 = parallel_ambiguity_analysis(ambiguity_sigmas, baseline_sigma, ambiguity_lambdas, baseline_lambda, 0.5, n_timesteps_short, g_map, n_iter=n_iter)

# Long term
ambiguity_long_sd1  = parallel_ambiguity_analysis(ambiguity_sigmas, baseline_sigma, ambiguity_lambdas, baseline_lambda, 1.0, n_timesteps_long, g_map, n_iter=n_iter)
ambiguity_long_sd2  = parallel_ambiguity_analysis(ambiguity_sigmas, baseline_sigma, ambiguity_lambdas, baseline_lambda, 2.0, n_timesteps_long, g_map, n_iter=n_iter)
ambiguity_long_sd05 = parallel_ambiguity_analysis(ambiguity_sigmas, baseline_sigma, ambiguity_lambdas, baseline_lambda, 0.5, n_timesteps_long, g_map, n_iter=n_iter)

# Save ambiguity results
CSV.write("data/ambiguity_results.csv", ambiguity_to_dataframe(vcat(
    ambiguity_short_sd1, ambiguity_short_sd2, ambiguity_short_sd05,
    ambiguity_long_sd1, ambiguity_long_sd2, ambiguity_long_sd05
)))
println("ambiguity_results.csv saved.")


#########################################
### SIGMA LAMBDA INTERACTION ANALYSIS ###
#########################################
interaction_sigmas = -4.2:.2:3.6
interaction_lambdas = -3.1:.2:4.7

# Unambiguous
sigma_lambda_short_DBM6 = parallel_sigma_lambda_analysis(interaction_sigmas, interaction_lambdas, 6.0, n_timesteps_short, g_map, n_iter=n_iter)
sigma_lambda_long_DBM6  = parallel_sigma_lambda_analysis(interaction_sigmas, interaction_lambdas, 6.0, n_timesteps_long, g_map, n_iter=n_iter)

# Semi-ambiguous
sigma_lambda_short_DBM45 = parallel_sigma_lambda_analysis(interaction_sigmas, interaction_lambdas, 4.5, n_timesteps_short, g_map, n_iter=n_iter)
sigma_lambda_long_DBM45  = parallel_sigma_lambda_analysis(interaction_sigmas, interaction_lambdas, 4.5, n_timesteps_long, g_map, n_iter=n_iter)

# Ambiguous
sigma_lambda_short_DBM3 = parallel_sigma_lambda_analysis(interaction_sigmas, interaction_lambdas, 3.0, n_timesteps_short, g_map, n_iter=n_iter)
sigma_lambda_long_DBM3  = parallel_sigma_lambda_analysis(interaction_sigmas, interaction_lambdas, 3.0, n_timesteps_long, g_map, n_iter=n_iter)

# Save interaction results
CSV.write("data/sigma_lambda_results.csv", sigma_lambda_to_dataframe(vcat(
    sigma_lambda_short_DBM6, sigma_lambda_long_DBM6,
    sigma_lambda_short_DBM45, sigma_lambda_long_DBM45,
    sigma_lambda_short_DBM3, sigma_lambda_long_DBM3
)))
println("sigma_lambda_results.csv saved.")


###########################
### TRAJECTORY ANALYSIS ###
###########################
trajectory_sigmas = 0.0:-0.8:-4
trajectory_lambdas = 0.7:-0.6:-2.3

# Unambiguous
trajectory_DBM6_span25  = parallel_trajectory_analysis(trajectory_sigmas, baseline_sigma, trajectory_lambdas, baseline_lambda, 6.0, 25, n_iter=n_iter)
trajectory_DBM6_span100 = parallel_trajectory_analysis(trajectory_sigmas, baseline_sigma, trajectory_lambdas, baseline_lambda, 6.0, 100, n_iter=n_iter)
trajectory_DBM6_span400 = parallel_trajectory_analysis(trajectory_sigmas, baseline_sigma, trajectory_lambdas, baseline_lambda, 6.0, 400, n_iter=n_iter)

# Ambiguous
trajectory_DBM3_span25  = parallel_trajectory_analysis(trajectory_sigmas, baseline_sigma, trajectory_lambdas, baseline_lambda, 3.0, 25, n_iter=n_iter)
trajectory_DBM3_span100 = parallel_trajectory_analysis(trajectory_sigmas, baseline_sigma, trajectory_lambdas, baseline_lambda, 3.0, 100, n_iter=n_iter)
trajectory_DBM3_span400 = parallel_trajectory_analysis(trajectory_sigmas, baseline_sigma, trajectory_lambdas, baseline_lambda, 3.0, 400, n_iter=n_iter)

# Save trajectory results
CSV.write("data/trajectory_results.csv", trajectory_to_dataframe(vcat(
    trajectory_DBM6_span25, trajectory_DBM6_span100, trajectory_DBM6_span400,
    trajectory_DBM3_span25, trajectory_DBM3_span100, trajectory_DBM3_span400
)))
println("trajectory_results.csv saved.")


#####################################
### PLOTTING TO CHECK THE RESULTS ###
#####################################
amb_df = CSV.read("data/ambiguity_results.csv", DataFrame)
sl_df  = CSV.read("data/sigma_lambda_results.csv", DataFrame)
tj_df  = CSV.read("data/trajectory_results.csv", DataFrame)

# Set metric here:
metric = :ncomp_mean  # :ari_mean, :mw2_mean, :ncomp_mean (also _min, _max, _q25, _median, _q75)

# SIGMA - short term
plot_ambiguity_analysis(filter(r -> r.lambda == 0.7 && r.n_timesteps == n_timesteps_short && r.sd == 1.0, amb_df), :sigma, "σ", metric=metric)
plot_ambiguity_analysis(filter(r -> r.lambda == 0.7 && r.n_timesteps == n_timesteps_short && r.sd == 2.0, amb_df), :sigma, "σ", metric=metric)
plot_ambiguity_analysis(filter(r -> r.lambda == 0.7 && r.n_timesteps == n_timesteps_short && r.sd == 0.5, amb_df), :sigma, "σ", metric=metric)

# SIGMA - long term
plot_ambiguity_analysis(filter(r -> r.lambda == 0.7 && r.n_timesteps == n_timesteps_long && r.sd == 1.0, amb_df), :sigma, "σ", metric=metric)
plot_ambiguity_analysis(filter(r -> r.lambda == 0.7 && r.n_timesteps == n_timesteps_long && r.sd == 2.0, amb_df), :sigma, "σ", metric=metric)
plot_ambiguity_analysis(filter(r -> r.lambda == 0.7 && r.n_timesteps == n_timesteps_long && r.sd == 0.5, amb_df), :sigma, "σ", metric=metric)

# LAMBDA - short term
plot_ambiguity_analysis(filter(r -> r.sigma == 0.0 && r.n_timesteps == n_timesteps_short && r.sd == 1.0, amb_df), :lambda, "λ", metric=metric)
plot_ambiguity_analysis(filter(r -> r.sigma == 0.0 && r.n_timesteps == n_timesteps_short && r.sd == 2.0, amb_df), :lambda, "λ", metric=metric)
plot_ambiguity_analysis(filter(r -> r.sigma == 0.0 && r.n_timesteps == n_timesteps_short && r.sd == 0.5, amb_df), :lambda, "λ", metric=metric)

# LAMBDA - long term
plot_ambiguity_analysis(filter(r -> r.sigma == 0.0 && r.n_timesteps == n_timesteps_long && r.sd == 1.0, amb_df), :lambda, "λ", metric=metric)
plot_ambiguity_analysis(filter(r -> r.sigma == 0.0 && r.n_timesteps == n_timesteps_long && r.sd == 2.0, amb_df), :lambda, "λ", metric=metric)
plot_ambiguity_analysis(filter(r -> r.sigma == 0.0 && r.n_timesteps == n_timesteps_long && r.sd == 0.5, amb_df), :lambda, "λ", metric=metric)

# Unambiguous interaction
plot_sigma_lambda_analysis(filter(r -> r.distance_between_means == 6.0 && r.n_timesteps == n_timesteps_short, sl_df), metric)
plot_sigma_lambda_analysis(filter(r -> r.distance_between_means == 6.0 && r.n_timesteps == n_timesteps_long, sl_df), metric)

# Semi-ambiguous interaction
plot_sigma_lambda_analysis(filter(r -> r.distance_between_means == 4.5 && r.n_timesteps == n_timesteps_short, sl_df), metric)
plot_sigma_lambda_analysis(filter(r -> r.distance_between_means == 4.5 && r.n_timesteps == n_timesteps_long, sl_df), metric)

# Ambiguous interaction
plot_sigma_lambda_analysis(filter(r -> r.distance_between_means == 3.0 && r.n_timesteps == n_timesteps_short, sl_df), metric)
plot_sigma_lambda_analysis(filter(r -> r.distance_between_means == 3.0 && r.n_timesteps == n_timesteps_long, sl_df), metric)

# Trajectories - sigma
ylims = startswith(String(metric), "ari") ? (.4,1) : :auto
plot_mean_trajectory(filter(r -> r.distance_between_means == 6.0 && r.timespan == 25  && r.lambda == 0.7, tj_df), :sigma, metric_mean=metric, ylims=ylims)
plot_mean_trajectory(filter(r -> r.distance_between_means == 6.0 && r.timespan == 100 && r.lambda == 0.7, tj_df), :sigma, metric_mean=metric, ylims=ylims)
plot_mean_trajectory(filter(r -> r.distance_between_means == 6.0 && r.timespan == 400 && r.lambda == 0.7, tj_df), :sigma, metric_mean=metric, ylims=ylims)

ylims = startswith(String(metric), "ari") ? (.4,.7) : :auto
plot_mean_trajectory(filter(r -> r.distance_between_means == 3.0 && r.timespan == 25  && r.lambda == 0.7, tj_df), :sigma, metric_mean=metric, ylims=ylims)
plot_mean_trajectory(filter(r -> r.distance_between_means == 3.0 && r.timespan == 100 && r.lambda == 0.7, tj_df), :sigma, metric_mean=metric, ylims=ylims)
plot_mean_trajectory(filter(r -> r.distance_between_means == 3.0 && r.timespan == 400 && r.lambda == 0.7, tj_df), :sigma, metric_mean=metric, ylims=ylims)

# Trajectories - lambda
ylims = startswith(String(metric), "ari") ? (.4,1) : :auto
plot_mean_trajectory(filter(r -> r.distance_between_means == 6.0 && r.timespan == 25  && r.sigma == 0.0, tj_df), :lambda, metric_mean=metric, ylims=ylims)
plot_mean_trajectory(filter(r -> r.distance_between_means == 6.0 && r.timespan == 100 && r.sigma == 0.0, tj_df), :lambda, metric_mean=metric, ylims=ylims)
plot_mean_trajectory(filter(r -> r.distance_between_means == 6.0 && r.timespan == 400 && r.sigma == 0.0, tj_df), :lambda, metric_mean=metric, ylims=ylims)

ylims = startswith(String(metric), "ari") ? (.4,.7) : :auto
plot_mean_trajectory(filter(r -> r.distance_between_means == 3.0 && r.timespan == 25  && r.sigma == 0.0, tj_df), :lambda, metric_mean=metric, ylims=ylims)
plot_mean_trajectory(filter(r -> r.distance_between_means == 3.0 && r.timespan == 100 && r.sigma == 0.0, tj_df), :lambda, metric_mean=metric, ylims=ylims)
plot_mean_trajectory(filter(r -> r.distance_between_means == 3.0 && r.timespan == 400 && r.sigma == 0.0, tj_df), :lambda, metric_mean=metric, ylims=ylims)
