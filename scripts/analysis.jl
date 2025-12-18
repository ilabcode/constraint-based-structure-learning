using JLD2

include("component.jl")  
include("mixture.jl")
include("simulation_functions.jl")  
include("analysis_functions.jl")    
include("plot_functions.jl")


##########################
### AMBIGUITY ANALYSIS ###
##########################
ambiguity_sigmas = -2.2:.2:2
ambiguity_lambdas = -2.5:.3:4 
n_timesteps_short = 100
n_timesteps_long = 10000

perfect_agent_ambiguity_results = ambiguity_analysis(perfect_agent = true)
jldsave("julia_data/perfect_agent_ambiguity_results.jld2"; perfect_agent_ambiguity_results)
println("perfect_agent_ambiguity_results have been saved.")

VERY_long_term_ambiguity_results = ambiguity_analysis(0.0,0.7,1.0,100000)
jldsave("julia_data/VERY_long_term_ambiguity_results.jld2"; VERY_long_term_ambiguity_results)
println("VERY_long_term_ambiguity_results have been saved.")

# SIGMA - short term 
sigma_ambiguity_short_sd1 = parallel_ambiguity_sigma(ambiguity_sigmas, 1.0, n_timesteps_short)
sigma_ambiguity_short_sd2 = parallel_ambiguity_sigma(ambiguity_sigmas, 2.0, n_timesteps_short)
sigma_ambiguity_short_sd05 = parallel_ambiguity_sigma(ambiguity_sigmas, 0.5, n_timesteps_short)
jldsave("julia_data/sigma_ambiguity_short.jld2"; sigma_ambiguity_short_sd1, sigma_ambiguity_short_sd2, sigma_ambiguity_short_sd05)
println("sigma_ambiguity_short have been saved.")


# SIGMA - long term
sigma_ambiguity_long_sd1 = parallel_ambiguity_sigma(ambiguity_sigmas, 1.0, n_timesteps_long)
sigma_ambiguity_long_sd2 = parallel_ambiguity_sigma(ambiguity_sigmas, 2.0, n_timesteps_long)
sigma_ambiguity_long_sd05 = parallel_ambiguity_sigma(ambiguity_sigmas, 0.5, n_timesteps_long)
jldsave("julia_data/sigma_ambiguity_long.jld2"; sigma_ambiguity_long_sd1, sigma_ambiguity_long_sd2, sigma_ambiguity_long_sd05)
println("sigma_ambiguity_long have been saved.")

# LAMBDA - short term
lambda_ambiguity_short_sd1 = parallel_ambiguity_lambda(ambiguity_lambdas, 1.0, n_timesteps_short)
lambda_ambiguity_short_sd2 = parallel_ambiguity_lambda(ambiguity_lambdas, 2.0, n_timesteps_short)
lambda_ambiguity_short_sd05 = parallel_ambiguity_lambda(ambiguity_lambdas, 0.5, n_timesteps_short)
jldsave("julia_data/lambda_ambiguity_short.jld2"; lambda_ambiguity_short_sd1, lambda_ambiguity_short_sd2, lambda_ambiguity_short_sd05)
println("lambda_ambiguity_short have been saved.")

# LAMBDA - long term
lambda_ambiguity_long_sd1 = parallel_ambiguity_lambda(ambiguity_lambdas, 1.0, n_timesteps_long)
lambda_ambiguity_long_sd2 = parallel_ambiguity_lambda(ambiguity_lambdas, 2.0, n_timesteps_long)
lambda_ambiguity_long_sd05 = parallel_ambiguity_lambda(ambiguity_lambdas, 0.5, n_timesteps_long)
jldsave("julia_data/lambda_ambiguity_long.jld2"; lambda_ambiguity_long_sd1, lambda_ambiguity_long_sd2, lambda_ambiguity_long_sd05)
println("lambda_ambiguity_long have been saved.")


#########################################
### SIGMA LAMBDA INTERACTION ANALYSIS ###
#########################################
interaction_sigmas = -4.2:.2:3.6
interaction_lambdas = -3.1:.2:4.7

# Unambiguous
sigma_lambda_short_DBM6 = parallel_sigma_lambda_analysis(interaction_sigmas, interaction_lambdas, 6.0, n_timesteps_short)
sigma_lambda_long_DBM6 = parallel_sigma_lambda_analysis(interaction_sigmas, interaction_lambdas, 6.0, n_timesteps_long)
jldsave("julia_data/sigma_lambda_DBM6.jld2"; sigma_lambda_short_DBM6, sigma_lambda_long_DBM6)
println("sigma_lambda_DBM6 have been saved.")

# Semi-ambiguous
sigma_lambda_short_DBM45 = parallel_sigma_lambda_analysis(interaction_sigmas, interaction_lambdas, 4.5, n_timesteps_short)
sigma_lambda_long_DBM45 = parallel_sigma_lambda_analysis(interaction_sigmas, interaction_lambdas, 4.5, n_timesteps_long)
jldsave("julia_data/sigma_lambda_DBM45.jld2"; sigma_lambda_short_DBM45, sigma_lambda_long_DBM45)
println("sigma_lambda_DBM45 have been saved.")

# Ambiguous
sigma_lambda_short_DBM3 = parallel_sigma_lambda_analysis(interaction_sigmas, interaction_lambdas, 3.0, n_timesteps_short)
sigma_lambda_long_DBM3 = parallel_sigma_lambda_analysis(interaction_sigmas, interaction_lambdas, 3.0, n_timesteps_long)
jldsave("julia_data/sigma_lambda_DBM3.jld2"; sigma_lambda_short_DBM3, sigma_lambda_long_DBM3)
println("sigma_lambda_DBM3 have been saved.")


###########################
### TRAJECTORY ANALYSIS ###
###########################
trajectory_sigmas = 0.0:-0.8:-4
trajectory_lambdas = 0.7:-0.6:-2.3

# SIGMA
# Unambiguous
sigma_trajectory_DBM6_span25 = parallel_trajectory_sigma(trajectory_sigmas, 6.0, 25)
sigma_trajectory_DBM6_span100 = parallel_trajectory_sigma(trajectory_sigmas, 6.0, 100)
sigma_trajectory_DBM6_span400 = parallel_trajectory_sigma(trajectory_sigmas, 6.0, 400)
jldsave("julia_data/sigma_trajectory_DBM6.jld2"; sigma_trajectory_DBM6_span25, sigma_trajectory_DBM6_span100, sigma_trajectory_DBM6_span400)
println("sigma_trajectory_DBM6 have been saved.")

# Ambiguous
sigma_trajectory_DBM3_span25 = parallel_trajectory_sigma(trajectory_sigmas, 3.0, 25)
sigma_trajectory_DBM3_span100 = parallel_trajectory_sigma(trajectory_sigmas, 3.0, 100)
sigma_trajectory_DBM3_span400 = parallel_trajectory_sigma(trajectory_sigmas, 3.0, 400)
jldsave("julia_data/sigma_trajectory_DBM3.jld2"; sigma_trajectory_DBM3_span25, sigma_trajectory_DBM3_span100, sigma_trajectory_DBM3_span400)
println("sigma_trajectory_DBM3 have been saved.")

# LAMBDA
# Unambiguous
lambda_trajectory_DBM6_span25 = parallel_trajectory_lambda(trajectory_lambdas, 6.0, 25)
lambda_trajectory_DBM6_span100 = parallel_trajectory_lambda(trajectory_lambdas, 6.0, 100)
lambda_trajectory_DBM6_span400 = parallel_trajectory_lambda(trajectory_lambdas, 6.0, 400)
jldsave("julia_data/lambda_trajectory_DBM6.jld2"; lambda_trajectory_DBM6_span25, lambda_trajectory_DBM6_span100, lambda_trajectory_DBM6_span400)
println("lambda_trajectory_DBM6 have been saved.")

# Ambiguous
lambda_trajectory_DBM3_span25 = parallel_trajectory_lambda(trajectory_lambdas, 3.0, 25)
lambda_trajectory_DBM3_span100 = parallel_trajectory_lambda(trajectory_lambdas, 3.0, 100)
lambda_trajectory_DBM3_span400 = parallel_trajectory_lambda(trajectory_lambdas, 3.0, 400)
jldsave("julia_data/lambda_trajectory_DBM3.jld2"; lambda_trajectory_DBM3_span25, lambda_trajectory_DBM3_span100, lambda_trajectory_DBM3_span400)
println("lambda_trajectory_DBM3 have been saved.")


#########################
### SAVE MASTER FILE ####
#########################

# jldsave("julia_data/all_analysis_results.jld2";
#     ambiguity_sigmas, ambiguity_lambdas, n_timesteps_short, n_timesteps_long,
#     perfect_agent_ambiguity_results, VERY_long_term_ambiguity_results,

#     # SIGMA - short term 
#     sigma_ambiguity_short_sd1, sigma_ambiguity_short_sd2, sigma_ambiguity_short_sd05,

#     # SIGMA - long term
#     sigma_ambiguity_long_sd1, sigma_ambiguity_long_sd2, sigma_ambiguity_long_sd05,

#     # LAMBDA - short term
#     lambda_ambiguity_short_sd1, lambda_ambiguity_short_sd2, lambda_ambiguity_short_sd05,

#     # LAMBDA - long term
#     lambda_ambiguity_long_sd1, lambda_ambiguity_long_sd2, lambda_ambiguity_long_sd05,

#     interaction_sigmas, interaction_lambdas,

#     # SIGMA LAMBDA - Unambiguous
#     sigma_lambda_short_DBM6, sigma_lambda_long_DBM6,

#     # SIGMA LAMBDA - Semi-ambiguous
#     sigma_lambda_short_DBM45, sigma_lambda_long_DBM45,

#     # SIGMA LAMBDA - Ambiguous
#     sigma_lambda_short_DBM3, sigma_lambda_long_DBM3,

#     trajectory_sigmas, trajectory_lambdas,

#     # SIGMA Trajectories
#     # Unambiguous
#     sigma_trajectory_DBM6_span25, sigma_trajectory_DBM6_span100, sigma_trajectory_DBM6_span400,

#     # Ambiguous
#     sigma_trajectory_DBM3_span25, sigma_trajectory_DBM3_span100, sigma_trajectory_DBM3_span400,

#     # LAMBDA Trajectories
#     # Unambiguous
#     lambda_trajectory_DBM6_span25, lambda_trajectory_DBM6_span100, lambda_trajectory_DBM6_span400,

#     # Ambiguous
#     lambda_trajectory_DBM3_span25, lambda_trajectory_DBM3_span100, lambda_trajectory_DBM3_span400
# )
# println("Aaaand... All done :)")


#####################################
### PLOTTING TO CHECK THE RESULTS ###
#####################################

# SIGMA - short term 
plot_ambiguity_analysis(sigma_ambiguity_short_sd1, ambiguity_sigmas, "σ")
plot_ambiguity_analysis(sigma_ambiguity_short_sd2, ambiguity_sigmas, "σ")
plot_ambiguity_analysis(sigma_ambiguity_short_sd05, ambiguity_sigmas, "σ")

# SIGMA - long term
plot_ambiguity_analysis(sigma_ambiguity_long_sd1, ambiguity_sigmas, "σ")
plot_ambiguity_analysis(sigma_ambiguity_long_sd2, ambiguity_sigmas, "σ")
plot_ambiguity_analysis(sigma_ambiguity_long_sd05, ambiguity_sigmas, "σ")

# LAMBDA - short term
plot_ambiguity_analysis(lambda_ambiguity_short_sd1, ambiguity_lambdas, "λ")
plot_ambiguity_analysis(lambda_ambiguity_short_sd2, ambiguity_lambdas, "λ")
plot_ambiguity_analysis(lambda_ambiguity_short_sd05, ambiguity_lambdas, "λ")

# LAMBDA - long term
plot_ambiguity_analysis(lambda_ambiguity_long_sd1, ambiguity_lambdas, "λ")
plot_ambiguity_analysis(lambda_ambiguity_long_sd2, ambiguity_lambdas, "λ")
plot_ambiguity_analysis(lambda_ambiguity_long_sd05, ambiguity_lambdas, "λ")

# Unambiguous interaction
plot_sigma_lambda_analysis(sigma_lambda_short_DBM6, interaction_sigmas, interaction_lambdas)
plot_sigma_lambda_analysis(sigma_lambda_long_DBM6, interaction_sigmas, interaction_lambdas)

# Semi-ambiguous interaction
plot_sigma_lambda_analysis(sigma_lambda_short_DBM45, interaction_sigmas, interaction_lambdas)
plot_sigma_lambda_analysis(sigma_lambda_long_DBM45, interaction_sigmas, interaction_lambdas)

# Ambiguous interaction
plot_sigma_lambda_analysis(sigma_lambda_short_DBM3, interaction_sigmas, interaction_lambdas)
plot_sigma_lambda_analysis(sigma_lambda_long_DBM3, interaction_sigmas, interaction_lambdas)


# Trajectories
plot_mean_trajectory(sigma_trajectory_DBM6_span25, "sigma", (.5,1))
plot_mean_trajectory(sigma_trajectory_DBM6_span100, "sigma", (.5,1))
plot_mean_trajectory(sigma_trajectory_DBM6_span400, "sigma", (.5,1))

plot_mean_trajectory(sigma_trajectory_DBM3_span25, "sigma", (.4,.7))
plot_mean_trajectory(sigma_trajectory_DBM3_span100, "sigma", (.4,.7))
plot_mean_trajectory(sigma_trajectory_DBM3_span400, "sigma", (.4,.7))

plot_mean_trajectory(lambda_trajectory_DBM6_span25, "lambda", (.4,1))
plot_mean_trajectory(lambda_trajectory_DBM6_span100, "lambda", (.4,1))
plot_mean_trajectory(lambda_trajectory_DBM6_span400, "lambda", (.4,1))

plot_mean_trajectory(lambda_trajectory_DBM3_span25, "lambda", (.4,.7))
plot_mean_trajectory(lambda_trajectory_DBM3_span100, "lambda", (.4,.7))
plot_mean_trajectory(lambda_trajectory_DBM3_span400, "lambda", (.4,.7))
