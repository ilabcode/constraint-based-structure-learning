using DataFrames, CSV
include("../plot_functions.jl")

data_dir = joinpath(@__DIR__, "../../data")
plots_dir = joinpath(@__DIR__, "../../plots/robustness analysis")

for gmap in ["logistic", "exp_decay"]

    save_plot(p, name) = (savefig(p, joinpath(plots_dir, "$(name)_$(gmap).png"));
                          savefig(p, joinpath(plots_dir, "$(name)_$(gmap).svg")))

    # Load analysis results
    amb_df     = CSV.read(joinpath(data_dir, "ambiguity_results_$(gmap).csv"), DataFrame)
    perfect_df = CSV.read(joinpath(data_dir, "perfect_agent_results_$(gmap).csv"), DataFrame)
    conv_df    = CSV.read(joinpath(data_dir, "convergence_baseline_results_$(gmap).csv"), DataFrame)
    sl_df      = CSV.read(joinpath(data_dir, "sigma_lambda_results_$(gmap).csv"), DataFrame)
    tj_df      = CSV.read(joinpath(data_dir, "trajectory_results_$(gmap).csv"), DataFrame)

    max_DBM6  = conv_df[conv_df.DBM .== 6.0, :ari_mean][1]
    max_DBM45 = conv_df[conv_df.DBM .== 4.5, :ari_mean][1]
    max_DBM3  = conv_df[conv_df.DBM .== 3.0, :ari_mean][1]

    conv_ari_DBM6   = conv_df[conv_df.DBM .== 6.0, :ari_mean][1]
    conv_ari_DBM3   = conv_df[conv_df.DBM .== 3.0, :ari_mean][1]
    conv_mw2_DBM6   = conv_df[conv_df.DBM .== 6.0, :mw2_mean][1]
    conv_mw2_DBM3   = conv_df[conv_df.DBM .== 3.0, :mw2_mean][1]
    conv_ncomp_DBM6 = conv_df[conv_df.DBM .== 6.0, :ncomp_mean][1]
    conv_ncomp_DBM3 = conv_df[conv_df.DBM .== 3.0, :ncomp_mean][1]

    # Ambiguity
    ambiguity_layout = @layout [ a{0.4w} b{0.4w} c{0.2w} ]

    sigma_ambiguity_plot_ari = plot(
        plot_ambiguity_analysis(filter(r -> r.lambda == 0.7 && r.n_timesteps == 100 && r.sd == 1.0, amb_df), :sigma, "σ",
            metric=:ari_mean, line_x=12, title="Short-term: 100 timesteps"),
        plot_ambiguity_analysis(filter(r -> r.lambda == 0.7 && r.n_timesteps == 10000 && r.sd == 1.0, amb_df), :sigma, "σ",
            metric=:ari_mean, line_x=12, title="Long-term: 10k timesteps"),
        plot_ambiguity_single(perfect_df, "Perfect agent", :ari_mean),

        layout = ambiguity_layout,
        size = (1000, 280),
        titlelocation = :left
    )

    mw2_clims = (0, 4)

    mw2_layout = @layout [ a{0.39w} b{0.5w} c{0.15w}  ]

    sigma_ambiguity_plot_mw2 = plot(
        plot_ambiguity_analysis(filter(r -> r.lambda == 0.7 && r.n_timesteps == 100 && r.sd == 1.0, amb_df), :sigma, "σ",
            metric=:mw2_mean, line_x=12, title="Short-term: 100 timesteps", reverse_cmap=true, clims=mw2_clims, text_color=:black),
        plot_ambiguity_analysis(filter(r -> r.lambda == 0.7 && r.n_timesteps == 10000 && r.sd == 1.0, amb_df), :sigma, "σ",
            metric=:mw2_mean, line_x=12, title="Long-term: 10k timesteps", reverse_cmap=true, clims=mw2_clims, colorbar=true, text_color=:black),
        layout = mw2_layout,
        size = (1000, 280),
        titlelocation = :left
    )

    lambda_ambiguity_plot_ari = plot(
        plot_ambiguity_analysis(filter(r -> r.sigma == 0.0 && r.n_timesteps == 100 && r.sd == 1.0, amb_df), :lambda, "λ",
            metric=:ari_mean, title="Short-term: 100 timesteps"),
        plot_ambiguity_analysis(filter(r -> r.sigma == 0.0 && r.n_timesteps == 10000 && r.sd == 1.0, amb_df), :lambda, "λ",
            metric=:ari_mean, title="Long-term: 10k timesteps"),
        plot_ambiguity_single(perfect_df, "Perfect agent", :ari_mean),

        layout = ambiguity_layout,
        size = (1000, 280),
        titlelocation = :left
    )

    lambda_ambiguity_plot_mw2 = plot(
        plot_ambiguity_analysis(filter(r -> r.sigma == 0.0 && r.n_timesteps == 100 && r.sd == 1.0, amb_df), :lambda, "λ",
            metric=:mw2_mean, title="Short-term: 100 timesteps", reverse_cmap=true, clims=mw2_clims, text_color=:black),
        plot_ambiguity_analysis(filter(r -> r.sigma == 0.0 && r.n_timesteps == 10000 && r.sd == 1.0, amb_df), :lambda, "λ",
            metric=:mw2_mean, title="Long-term: 10k timesteps", reverse_cmap=true, clims=mw2_clims, colorbar=true, text_color=:black),
        layout = mw2_layout,
        size = (1000, 280),
        titlelocation = :left
    )

    # Number of components (deviation from optimal = 3)
    amb_df.ncomp_dev_mean = abs.(amb_df.ncomp_mean .- 3.0)
    ncomp_clims = (0, 4)

    sigma_ambiguity_plot_ncomp = plot(
        plot_ambiguity_analysis(filter(r -> r.lambda == 0.7 && r.n_timesteps == 100 && r.sd == 1.0, amb_df), :sigma, "σ",
            metric=:ncomp_dev_mean, line_x=12, title="Short-term: 100 timesteps", reverse_cmap=true, clims=ncomp_clims, text_color=:black),
        plot_ambiguity_analysis(filter(r -> r.lambda == 0.7 && r.n_timesteps == 10000 && r.sd == 1.0, amb_df), :sigma, "σ",
            metric=:ncomp_dev_mean, line_x=12, title="Long-term: 10k timesteps", reverse_cmap=true, clims=ncomp_clims, colorbar=true, text_color=:black),
        layout = mw2_layout,
        size = (1000, 280),
        titlelocation = :left
    )

    lambda_ambiguity_plot_ncomp = plot(
        plot_ambiguity_analysis(filter(r -> r.sigma == 0.0 && r.n_timesteps == 100 && r.sd == 1.0, amb_df), :lambda, "λ",
            metric=:ncomp_dev_mean, title="Short-term: 100 timesteps", reverse_cmap=true, clims=ncomp_clims, text_color=:black),
        plot_ambiguity_analysis(filter(r -> r.sigma == 0.0 && r.n_timesteps == 10000 && r.sd == 1.0, amb_df), :lambda, "λ",
            metric=:ncomp_dev_mean, title="Long-term: 10k timesteps", reverse_cmap=true, clims=ncomp_clims, colorbar=true, text_color=:black),
        layout = mw2_layout,
        size = (1000, 280),
        titlelocation = :left
    )

    # Sigma lambda interaction
    interaction_layout = @layout [ a{0.425w} b{0.585w}
                                   c{0.425w} d{0.585w}
                                   e{0.425w} f{0.585w} ]

    # Per-row clims: max ARI in the long-term (n_timesteps=10000) data at each DBM
    lt_sl = filter(r -> r.n_timesteps == 10000, sl_df)
    row_max_DBM6  = maximum(filter(r -> r.distance_between_means == 6.0, lt_sl).ari_mean)
    row_max_DBM45 = maximum(filter(r -> r.distance_between_means == 4.5, lt_sl).ari_mean)
    row_max_DBM3  = maximum(filter(r -> r.distance_between_means == 3.0, lt_sl).ari_mean)

    interaction_plot = plot(
        plot_sigma_lambda_analysis(filter(r -> r.distance_between_means == 6.0 && r.n_timesteps == 100, sl_df), :ari_mean,
            clims=(0, row_max_DBM6), legend=false, xlab="", mark_baseline=true, mark_best=true),
        plot_sigma_lambda_analysis(filter(r -> r.distance_between_means == 6.0 && r.n_timesteps == 10000, sl_df), :ari_mean,
            clims=(0, row_max_DBM6), ylab="", xlab="", mark_baseline=true, mark_best=true),
        plot_sigma_lambda_analysis(filter(r -> r.distance_between_means == 4.5 && r.n_timesteps == 100, sl_df), :ari_mean,
            clims=(0, row_max_DBM45), legend=false, xlab="", mark_baseline=true, mark_best=true),
        plot_sigma_lambda_analysis(filter(r -> r.distance_between_means == 4.5 && r.n_timesteps == 10000, sl_df), :ari_mean,
            clims=(0, row_max_DBM45), ylab="", xlab="", mark_baseline=true, mark_best=true),
        plot_sigma_lambda_analysis(filter(r -> r.distance_between_means == 3.0 && r.n_timesteps == 100, sl_df), :ari_mean,
            clims=(0, row_max_DBM3), legend=false, mark_baseline=true, mark_best=true),
        plot_sigma_lambda_analysis(filter(r -> r.distance_between_means == 3.0 && r.n_timesteps == 10000, sl_df), :ari_mean,
            clims=(0, row_max_DBM3), ylab="", mark_baseline=true, mark_best=true),
        layout = interaction_layout,
        size = (800, 900)
    )

    # Sigma lambda interaction — MW2
    interaction_plot_mw2 = plot(
        plot_sigma_lambda_analysis(filter(r -> r.distance_between_means == 6.0 && r.n_timesteps == 100, sl_df), :mw2_mean,
            clims=:auto, legend=false, xlab="", reverse_cmap=true),
        plot_sigma_lambda_analysis(filter(r -> r.distance_between_means == 6.0 && r.n_timesteps == 10000, sl_df), :mw2_mean,
            clims=:auto, ylab="", xlab="", reverse_cmap=true),
        plot_sigma_lambda_analysis(filter(r -> r.distance_between_means == 4.5 && r.n_timesteps == 100, sl_df), :mw2_mean,
            clims=:auto, legend=false, xlab="", reverse_cmap=true),
        plot_sigma_lambda_analysis(filter(r -> r.distance_between_means == 4.5 && r.n_timesteps == 10000, sl_df), :mw2_mean,
            clims=:auto, ylab="", xlab="", reverse_cmap=true),
        plot_sigma_lambda_analysis(filter(r -> r.distance_between_means == 3.0 && r.n_timesteps == 100, sl_df), :mw2_mean,
            clims=:auto, legend=false, reverse_cmap=true),
        plot_sigma_lambda_analysis(filter(r -> r.distance_between_means == 3.0 && r.n_timesteps == 10000, sl_df), :mw2_mean,
            clims=:auto, ylab="", reverse_cmap=true),
        layout = interaction_layout,
        size = (800, 900)
    )

    # Sigma lambda interaction — N components (deviation from 3)
    sl_df.ncomp_dev_mean = abs.(sl_df.ncomp_mean .- 3.0)

    interaction_plot_ncomp = plot(
        plot_sigma_lambda_analysis(filter(r -> r.distance_between_means == 6.0 && r.n_timesteps == 100, sl_df), :ncomp_dev_mean,
            clims=:auto, legend=false, xlab="", reverse_cmap=true),
        plot_sigma_lambda_analysis(filter(r -> r.distance_between_means == 6.0 && r.n_timesteps == 10000, sl_df), :ncomp_dev_mean,
            clims=:auto, ylab="", xlab="", reverse_cmap=true),
        plot_sigma_lambda_analysis(filter(r -> r.distance_between_means == 4.5 && r.n_timesteps == 100, sl_df), :ncomp_dev_mean,
            clims=:auto, legend=false, xlab="", reverse_cmap=true),
        plot_sigma_lambda_analysis(filter(r -> r.distance_between_means == 4.5 && r.n_timesteps == 10000, sl_df), :ncomp_dev_mean,
            clims=:auto, ylab="", xlab="", reverse_cmap=true),
        plot_sigma_lambda_analysis(filter(r -> r.distance_between_means == 3.0 && r.n_timesteps == 100, sl_df), :ncomp_dev_mean,
            clims=:auto, legend=false, reverse_cmap=true),
        plot_sigma_lambda_analysis(filter(r -> r.distance_between_means == 3.0 && r.n_timesteps == 10000, sl_df), :ncomp_dev_mean,
            clims=:auto, ylab="", reverse_cmap=true),
        layout = interaction_layout,
        size = (800, 900)
    )

    # Trajectories — 3×3 grids: rows = metrics (ARI, MW2, Ncomp), cols = spans (25, 100, 400)
    sigma_subset  = [0.0, -1.6, -2.4, -4.0]
    lambda_subset = [0.7, 0.1, -0.5, -1.7]

    # Sigma — DBM 6
    sigma_trajectory_DBM6 = plot(
        plot_mean_trajectory(filter(r -> r.distance_between_means == 6.0 && r.timespan == 25  && r.lambda == 0.7 && r.sigma in sigma_subset, tj_df), :sigma, metric_mean=:ari_mean, ylims=(.4,1), perfect_level=conv_ari_DBM6, legend=:bottomright),
        plot_mean_trajectory(filter(r -> r.distance_between_means == 6.0 && r.timespan == 100 && r.lambda == 0.7 && r.sigma in sigma_subset, tj_df), :sigma, metric_mean=:ari_mean, ylims=(.4,1), perfect_level=conv_ari_DBM6),
        plot_mean_trajectory(filter(r -> r.distance_between_means == 6.0 && r.timespan == 400 && r.lambda == 0.7 && r.sigma in sigma_subset, tj_df), :sigma, metric_mean=:ari_mean, ylims=(.4,1), perfect_level=conv_ari_DBM6),

        plot_mean_trajectory(filter(r -> r.distance_between_means == 6.0 && r.timespan == 25  && r.lambda == 0.7 && r.sigma in sigma_subset, tj_df), :sigma, metric_mean=:mw2_mean, ylims=(0.9, 2.2), perfect_level=conv_mw2_DBM6),
        plot_mean_trajectory(filter(r -> r.distance_between_means == 6.0 && r.timespan == 100 && r.lambda == 0.7 && r.sigma in sigma_subset, tj_df), :sigma, metric_mean=:mw2_mean, ylims=(0.9, 2.2), perfect_level=conv_mw2_DBM6),
        plot_mean_trajectory(filter(r -> r.distance_between_means == 6.0 && r.timespan == 400 && r.lambda == 0.7 && r.sigma in sigma_subset, tj_df), :sigma, metric_mean=:mw2_mean, ylims=(0.9, 2.2), perfect_level=conv_mw2_DBM6),

        plot_mean_trajectory(filter(r -> r.distance_between_means == 6.0 && r.timespan == 25  && r.lambda == 0.7 && r.sigma in sigma_subset, tj_df), :sigma, metric_mean=:ncomp_mean, ylims=:auto, perfect_level=conv_ncomp_DBM6),
        plot_mean_trajectory(filter(r -> r.distance_between_means == 6.0 && r.timespan == 100 && r.lambda == 0.7 && r.sigma in sigma_subset, tj_df), :sigma, metric_mean=:ncomp_mean, ylims=:auto, perfect_level=conv_ncomp_DBM6),
        plot_mean_trajectory(filter(r -> r.distance_between_means == 6.0 && r.timespan == 400 && r.lambda == 0.7 && r.sigma in sigma_subset, tj_df), :sigma, metric_mean=:ncomp_mean, ylims=:auto, perfect_level=conv_ncomp_DBM6),
        layout = (3, 3),
        size = (1000, 600)
    )

    # Sigma — DBM 3
    sigma_trajectory_DBM3 = plot(
        plot_mean_trajectory(filter(r -> r.distance_between_means == 3.0 && r.timespan == 25  && r.lambda == 0.7 && r.sigma in sigma_subset, tj_df), :sigma, metric_mean=:ari_mean, ylims=(.35,.8), perfect_level=conv_ari_DBM3, legend=:bottomright),
        plot_mean_trajectory(filter(r -> r.distance_between_means == 3.0 && r.timespan == 100 && r.lambda == 0.7 && r.sigma in sigma_subset, tj_df), :sigma, metric_mean=:ari_mean, ylims=(.35,.8), perfect_level=conv_ari_DBM3),
        plot_mean_trajectory(filter(r -> r.distance_between_means == 3.0 && r.timespan == 400 && r.lambda == 0.7 && r.sigma in sigma_subset, tj_df), :sigma, metric_mean=:ari_mean, ylims=(.35,.8), perfect_level=conv_ari_DBM3),

        plot_mean_trajectory(filter(r -> r.distance_between_means == 3.0 && r.timespan == 25  && r.lambda == 0.7 && r.sigma in sigma_subset, tj_df), :sigma, metric_mean=:mw2_mean, ylims=(0.9, 2.2), perfect_level=conv_mw2_DBM3),
        plot_mean_trajectory(filter(r -> r.distance_between_means == 3.0 && r.timespan == 100 && r.lambda == 0.7 && r.sigma in sigma_subset, tj_df), :sigma, metric_mean=:mw2_mean, ylims=(0.9, 2.2), perfect_level=conv_mw2_DBM3),
        plot_mean_trajectory(filter(r -> r.distance_between_means == 3.0 && r.timespan == 400 && r.lambda == 0.7 && r.sigma in sigma_subset, tj_df), :sigma, metric_mean=:mw2_mean, ylims=(0.9, 2.2), perfect_level=conv_mw2_DBM3),

        plot_mean_trajectory(filter(r -> r.distance_between_means == 3.0 && r.timespan == 25  && r.lambda == 0.7 && r.sigma in sigma_subset, tj_df), :sigma, metric_mean=:ncomp_mean, ylims=:auto, perfect_level=conv_ncomp_DBM3),
        plot_mean_trajectory(filter(r -> r.distance_between_means == 3.0 && r.timespan == 100 && r.lambda == 0.7 && r.sigma in sigma_subset, tj_df), :sigma, metric_mean=:ncomp_mean, ylims=:auto, perfect_level=conv_ncomp_DBM3),
        plot_mean_trajectory(filter(r -> r.distance_between_means == 3.0 && r.timespan == 400 && r.lambda == 0.7 && r.sigma in sigma_subset, tj_df), :sigma, metric_mean=:ncomp_mean, ylims=:auto, perfect_level=conv_ncomp_DBM3),
        layout = (3, 3),
        size = (1000, 600)
    )

    # Lambda — DBM 6
    lambda_trajectory_DBM6 = plot(
        plot_mean_trajectory(filter(r -> r.distance_between_means == 6.0 && r.timespan == 25  && r.sigma == 0.0 && r.lambda in lambda_subset, tj_df), :lambda, metric_mean=:ari_mean, ylims=(.3,1), perfect_level=conv_ari_DBM6, legend=:bottomright),
        plot_mean_trajectory(filter(r -> r.distance_between_means == 6.0 && r.timespan == 100 && r.sigma == 0.0 && r.lambda in lambda_subset, tj_df), :lambda, metric_mean=:ari_mean, ylims=(.3,1), perfect_level=conv_ari_DBM6),
        plot_mean_trajectory(filter(r -> r.distance_between_means == 6.0 && r.timespan == 400 && r.sigma == 0.0 && r.lambda in lambda_subset, tj_df), :lambda, metric_mean=:ari_mean, ylims=(.3,1), perfect_level=conv_ari_DBM6),

        plot_mean_trajectory(filter(r -> r.distance_between_means == 6.0 && r.timespan == 25  && r.sigma == 0.0 && r.lambda in lambda_subset, tj_df), :lambda, metric_mean=:mw2_mean, ylims=(0.9, 2.2), perfect_level=conv_mw2_DBM6),
        plot_mean_trajectory(filter(r -> r.distance_between_means == 6.0 && r.timespan == 100 && r.sigma == 0.0 && r.lambda in lambda_subset, tj_df), :lambda, metric_mean=:mw2_mean, ylims=(0.9, 2.2), perfect_level=conv_mw2_DBM6),
        plot_mean_trajectory(filter(r -> r.distance_between_means == 6.0 && r.timespan == 400 && r.sigma == 0.0 && r.lambda in lambda_subset, tj_df), :lambda, metric_mean=:mw2_mean, ylims=(0.9, 2.2), perfect_level=conv_mw2_DBM6),

        plot_mean_trajectory(filter(r -> r.distance_between_means == 6.0 && r.timespan == 25  && r.sigma == 0.0 && r.lambda in lambda_subset, tj_df), :lambda, metric_mean=:ncomp_mean, ylims=:auto, perfect_level=conv_ncomp_DBM6),
        plot_mean_trajectory(filter(r -> r.distance_between_means == 6.0 && r.timespan == 100 && r.sigma == 0.0 && r.lambda in lambda_subset, tj_df), :lambda, metric_mean=:ncomp_mean, ylims=:auto, perfect_level=conv_ncomp_DBM6),
        plot_mean_trajectory(filter(r -> r.distance_between_means == 6.0 && r.timespan == 400 && r.sigma == 0.0 && r.lambda in lambda_subset, tj_df), :lambda, metric_mean=:ncomp_mean, ylims=:auto, perfect_level=conv_ncomp_DBM6),
        layout = (3, 3),
        size = (1000, 600)
    )

    # Lambda — DBM 3
    lambda_trajectory_DBM3 = plot(
        plot_mean_trajectory(filter(r -> r.distance_between_means == 3.0 && r.timespan == 25  && r.sigma == 0.0 && r.lambda in lambda_subset, tj_df), :lambda, metric_mean=:ari_mean, ylims=(.35,.8), perfect_level=conv_ari_DBM3, legend=:bottomright),
        plot_mean_trajectory(filter(r -> r.distance_between_means == 3.0 && r.timespan == 100 && r.sigma == 0.0 && r.lambda in lambda_subset, tj_df), :lambda, metric_mean=:ari_mean, ylims=(.35,.8), perfect_level=conv_ari_DBM3),
        plot_mean_trajectory(filter(r -> r.distance_between_means == 3.0 && r.timespan == 400 && r.sigma == 0.0 && r.lambda in lambda_subset, tj_df), :lambda, metric_mean=:ari_mean, ylims=(.35,.8), perfect_level=conv_ari_DBM3),

        plot_mean_trajectory(filter(r -> r.distance_between_means == 3.0 && r.timespan == 25  && r.sigma == 0.0 && r.lambda in lambda_subset, tj_df), :lambda, metric_mean=:mw2_mean, ylims=(0.9, 2.2), perfect_level=conv_mw2_DBM3),
        plot_mean_trajectory(filter(r -> r.distance_between_means == 3.0 && r.timespan == 100 && r.sigma == 0.0 && r.lambda in lambda_subset, tj_df), :lambda, metric_mean=:mw2_mean, ylims=(0.9, 2.2), perfect_level=conv_mw2_DBM3),
        plot_mean_trajectory(filter(r -> r.distance_between_means == 3.0 && r.timespan == 400 && r.sigma == 0.0 && r.lambda in lambda_subset, tj_df), :lambda, metric_mean=:mw2_mean, ylims=(0.9, 2.2), perfect_level=conv_mw2_DBM3),

        plot_mean_trajectory(filter(r -> r.distance_between_means == 3.0 && r.timespan == 25  && r.sigma == 0.0 && r.lambda in lambda_subset, tj_df), :lambda, metric_mean=:ncomp_mean, ylims=:auto, perfect_level=conv_ncomp_DBM3),
        plot_mean_trajectory(filter(r -> r.distance_between_means == 3.0 && r.timespan == 100 && r.sigma == 0.0 && r.lambda in lambda_subset, tj_df), :lambda, metric_mean=:ncomp_mean, ylims=:auto, perfect_level=conv_ncomp_DBM3),
        plot_mean_trajectory(filter(r -> r.distance_between_means == 3.0 && r.timespan == 400 && r.sigma == 0.0 && r.lambda in lambda_subset, tj_df), :lambda, metric_mean=:ncomp_mean, ylims=:auto, perfect_level=conv_ncomp_DBM3),
        layout = (3, 3),
        size = (1000, 600)
    )

    # Distribution comparison (baseline parameters: sigma=0, lambda=0.7)
    baseline_amb = filter(r -> r.sigma == 0.0 && r.lambda == 0.7 && r.sd == 1.0, amb_df)

    distribution_plot = plot(
        plot_distribution_comparison(baseline_amb, nothing, perfect_df, "ari",
            conv_df=conv_df, title="ARI", ylabel=" "),
        plot_distribution_comparison(baseline_amb, nothing, nothing, "mw2",
            conv_df=conv_df, title="MW₂", ylabel=" ", show_legend=false),
        plot_distribution_comparison(baseline_amb, nothing, nothing, "ncomp",
            conv_df=conv_df, title="Number of Components", ylabel=" ", show_legend=false),
        layout = (1, 3),
        size = (1500, 450)
    )

    # Distribution comparison for sigma
    low_sigma_amb = filter(r -> r.sigma == -1.0 && r.lambda == 0.7 && r.sd == 1.0, amb_df)
    high_sigma_amb = filter(r -> r.sigma == 2.0 && r.lambda == 0.7 && r.sd == 1.0, amb_df)

    sigma_distribution_plot = plot(
        plot_distribution_comparison(low_sigma_amb, baseline_amb, nothing, "ari",
            title="ARI", ylabel="σ = −1.0"),
        plot_distribution_comparison(low_sigma_amb, baseline_amb, nothing, "mw2",
            title="MW₂", ylabel=" ", show_legend=false),
        plot_distribution_comparison(low_sigma_amb, baseline_amb, nothing, "ncomp",
            title="Component Count", ylabel=" ", show_legend=false),
        plot_distribution_comparison(high_sigma_amb, baseline_amb, nothing, "ari",
            ylabel="σ = 2.0", show_legend=false),
        plot_distribution_comparison(high_sigma_amb, baseline_amb, nothing, "mw2",
            ylabel=" ", show_legend=false),
        plot_distribution_comparison(high_sigma_amb, baseline_amb, nothing, "ncomp",
            ylabel=" ", show_legend=false),
        layout = (2, 3),
        size = (1500, 850)
    )

    # Distribution comparison for lambda
    low_lambda_amb = filter(r -> r.sigma == 0.0 && r.lambda == -1.0 && r.sd == 1.0, amb_df)
    high_lambda_amb = filter(r -> r.sigma == 0.0 && r.lambda == 2.0 && r.sd == 1.0, amb_df)

    lambda_distribution_plot = plot(
        plot_distribution_comparison(low_lambda_amb, baseline_amb, nothing, "ari",
            title="ARI", ylabel="λ = −1.0"),
        plot_distribution_comparison(low_lambda_amb, baseline_amb, nothing, "mw2",
            title="MW₂", ylabel=" ", show_legend=false),
        plot_distribution_comparison(low_lambda_amb, baseline_amb, nothing, "ncomp",
            title="Component Count", ylabel=" ", show_legend=false),
        plot_distribution_comparison(high_lambda_amb, baseline_amb, nothing, "ari",
            ylabel="λ = 2.0", show_legend=false),
        plot_distribution_comparison(high_lambda_amb, baseline_amb, nothing, "mw2",
            ylabel=" ", show_legend=false),
        plot_distribution_comparison(high_lambda_amb, baseline_amb, nothing, "ncomp",
            ylabel=" ", show_legend=false),
        layout = (2, 3),
        size = (1500, 850)
    )

    # Save plots (PNG + SVG)
    # Ambiguity
    save_plot(sigma_ambiguity_plot_ari, "sigma_ambiguity_ari")
    save_plot(sigma_ambiguity_plot_mw2, "sigma_ambiguity_mw2")
    save_plot(sigma_ambiguity_plot_ncomp, "sigma_ambiguity_ncomp")
    save_plot(lambda_ambiguity_plot_ari, "lambda_ambiguity_ari")
    save_plot(lambda_ambiguity_plot_mw2, "lambda_ambiguity_mw2")
    save_plot(lambda_ambiguity_plot_ncomp, "lambda_ambiguity_ncomp")

    # Interaction
    save_plot(interaction_plot, "interaction_ari")
    save_plot(interaction_plot_mw2, "interaction_mw2")
    save_plot(interaction_plot_ncomp, "interaction_ncomp")

    # Trajectories
    save_plot(sigma_trajectory_DBM6, "sigma_trajectory_DBM6")
    save_plot(sigma_trajectory_DBM3, "sigma_trajectory_DBM3")
    save_plot(lambda_trajectory_DBM6, "lambda_trajectory_DBM6")
    save_plot(lambda_trajectory_DBM3, "lambda_trajectory_DBM3")

    # Distribution comparisons
    save_plot(distribution_plot, "distribution_baseline")
    save_plot(sigma_distribution_plot, "distribution_sigma")
    save_plot(lambda_distribution_plot, "distribution_lambda")

    println("Plots saved for $(gmap).")

end
