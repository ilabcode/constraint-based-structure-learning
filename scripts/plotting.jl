using JLD2
include("plot_functions.jl")

# Load analysis results
AR = load("julia_data/all_analysis_results.jld2")

max_DBM6 = AR["VERY_long_term_ambiguity_results"][12]
max_DBM45 = AR["VERY_long_term_ambiguity_results"][9]
max_DBM3 = AR["VERY_long_term_ambiguity_results"][6]

# Ambiguity
ambiguity_layout = @layout    [ a{0.4w}  b{0.4w} c{0.2w} ]
sigma_ambiguity_plot= plot(
    plot_ambiguity_analysis(AR["sigma_ambiguity_short_sd1"], AR["ambiguity_sigmas"], "σ", 12, "Short-term: 100 timesteps"),
    plot_ambiguity_analysis(AR["sigma_ambiguity_long_sd1"], AR["ambiguity_sigmas"], "σ", 12, "Long-term: 10k timesteps"),
    plot_ambiguity_single(AR["perfect_agent_ambiguity_results"], "Perfect agent"),

    layout = ambiguity_layout,
    size = (1000,280),
    titlelocation = :left         
)

lambda_ambiguity_plot= plot(
    plot_ambiguity_analysis(AR["lambda_ambiguity_short_sd1"], AR["ambiguity_lambdas"], "λ", nothing, "Short-term: 100 timesteps"),
    plot_ambiguity_analysis(AR["lambda_ambiguity_long_sd1"], AR["ambiguity_lambdas"], "λ", nothing, "Long-term: 10k timesteps"),
    plot_ambiguity_single(AR["perfect_agent_ambiguity_results"], "Perfect agent"),

    layout = ambiguity_layout,
    size = (1000,280),
    titlelocation = :left         
)


# Sigma lambda
interaction_layout = @layout    [ a{0.425w} b{0.585w}
                                c{0.425w} d{0.585w}
                                e{0.425w} f{0.585w}]
interaction_plot = plot(
    plot_sigma_lambda_analysis(AR["sigma_lambda_short_DBM6"], AR["interaction_sigmas"], AR["interaction_lambdas"], max_DBM6; legend=false, xlab=""),
    plot_sigma_lambda_analysis(AR["sigma_lambda_long_DBM6"], AR["interaction_sigmas"], AR["interaction_lambdas"], max_DBM6; ylab="", xlab=""),
    plot_sigma_lambda_analysis(AR["sigma_lambda_short_DBM45"], AR["interaction_sigmas"], AR["interaction_lambdas"], max_DBM45; legend=false, xlab=""),
    plot_sigma_lambda_analysis(AR["sigma_lambda_long_DBM45"], AR["interaction_sigmas"], AR["interaction_lambdas"], max_DBM45; ylab="", xlab=""),
    plot_sigma_lambda_analysis(AR["sigma_lambda_short_DBM3"], AR["interaction_sigmas"], AR["interaction_lambdas"], max_DBM3; legend=false),
    plot_sigma_lambda_analysis(AR["sigma_lambda_long_DBM3"], AR["interaction_sigmas"], AR["interaction_lambdas"], max_DBM3; ylab=""),
    layout=interaction_layout,
    size = (800,900)
)

# Trajectories

sigma_lines=[1,3,4,6]
lambda_lines=[1,2,3,5]
sigma_trajectory_plot = plot(
    plot_mean_trajectory(AR["sigma_trajectory_DBM6_span25"][sigma_lines], "sigma", (.4,1), max_DBM6, :bottomright),
    plot_mean_trajectory(AR["sigma_trajectory_DBM6_span100"][sigma_lines], "sigma", (.4,1), max_DBM6),
    plot_mean_trajectory(AR["sigma_trajectory_DBM6_span400"][sigma_lines], "sigma", (.4,1), max_DBM6),

    plot_mean_trajectory(AR["sigma_trajectory_DBM3_span25"][sigma_lines], "sigma", (.35,.8), max_DBM3),
    plot_mean_trajectory(AR["sigma_trajectory_DBM3_span100"][sigma_lines], "sigma", (.35,.8), max_DBM3),
    plot_mean_trajectory(AR["sigma_trajectory_DBM3_span400"][sigma_lines], "sigma", (.35,.8), max_DBM3),
    size=(1000,400)
)

lambda_trajectory_plot = plot(
    plot_mean_trajectory(AR["lambda_trajectory_DBM6_span25"][lambda_lines], "lambda", (.3,1), max_DBM6, :bottomright),
    plot_mean_trajectory(AR["lambda_trajectory_DBM6_span100"][lambda_lines], "lambda", (.3,1), max_DBM6),
    plot_mean_trajectory(AR["lambda_trajectory_DBM6_span400"][lambda_lines], "lambda", (.3,1), max_DBM6),

    plot_mean_trajectory(AR["lambda_trajectory_DBM3_span25"][lambda_lines], "lambda", (.35,.8), max_DBM3),
    plot_mean_trajectory(AR["lambda_trajectory_DBM3_span100"][lambda_lines], "lambda", (.35,.8), max_DBM3),
    plot_mean_trajectory(AR["lambda_trajectory_DBM3_span400"][lambda_lines], "lambda", (.35,.8), max_DBM3), 
    size=(1000,400)
)

# Save plots
# savefig(sigma_ambiguity_plot, "plots/sigma_ambiguity_plot.svg")
# savefig(lambda_ambiguity_plot, "plots/lambda_ambiguity_plot.svg")
# savefig(sigma_trajectory_plot, "plots/sigma_trajectory_plot.png")
# savefig(lambda_trajectory_plot, "plots/lambda_trajectory_plot.png")
# savefig(interaction_plot, "plots/interaction_plot.svg")

# Example Plots
example_plot(
    0.0, #Sigma
    0.7, #lambda
    10, #n_timesteps
    6, # Distance between means
    (200, 1000);
    ymax = .9
    )

appendix1 = many_simulations_plot(
    100, # N simulations
    0.0, #sigma
    0.7, #lambda
    get_simple_environment(6.0),
    (1200,1000), # size
    100, # n timesteps per simulation
    (10,10)
    )
savefig(appendix1, "plots/100_standard_DBM6_shortterm.png")

appendix2 = many_simulations_plot(
    100, # N simulations
    0.0, #sigma
    0.7, #lambda
    get_simple_environment(3.0),
    (1200,1000), # size
    10000, # n timesteps per simulation
    (10,10)
    )
savefig(appendix1, "plots/100_standard_DBM3_longterm.png")