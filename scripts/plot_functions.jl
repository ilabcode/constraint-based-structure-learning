using Plots; gr()
using Measures
using PerceptualColourMaps
using Statistics
using Smoothing
using DataFrames
include("component.jl")
include("mixture.jl")
include("simulation_functions.jl")

function plot_simulation_environments(thetas, n_timesteps; x_limits = [-10, 10], ymax = .5)
    plot_list = []
    x = x_limits[1]:.01:x_limits[2]
    for i in eachindex(thetas)
        plt = plot(ylim=[0,ymax], title = string(n_timesteps[i])*" timesteps")
        for theta in thetas[i] 
            gausian_data = [pdf(Normal(theta...),i) for i in x]
            gausian_data .*= 1/length(thetas[i] )
            plot!(x, zeros(length(x)), fillrange = gausian_data, c = "black", alpha = .2, label = nothing)
        end
        push!(plot_list, plt)
    end

    plot(plot_list...)
end

function plot_mixture(mixture; x_limits = [-10, 10], ymax = .5, title = "", environment = nothing, alpha = 1, colors = nothing, background=:white)
    
    x = x_limits[1]:.01:x_limits[2]
    plt = plot(ylim=[0,ymax], title = title, axis=false, background_color = background)

    if !isnothing(environment)
        for env in environment
            gausian_data = [pdf(Normal(env...),i) for i in x]
            plot!(x, zeros(length(x)), fillrange = gausian_data, c = "gray", alpha = 0.3, label = nothing, yaxis=false)
        end
    end

    for i in eachindex(mixture.components)
        component = mixture.components[i]
        gausian_data = [pdf(Normal(component.theta...),i) for i in x] .*= component.weight
        if isnothing(colors)
            color = mixture.IDs[i]
        else 
            if mixture.IDs[i] > length(colors)
                color = :grey
            else
                color = colors[mixture.IDs[i]]
            end
        end

        plot!(x, gausian_data, linewidth=2, c = color, label = nothing, grid = false, alpha = alpha)
    end

    return plt
end

function many_simulations_plot(n_sim, sigma, lambda, environment, size = (1000,700), n_timesteps=100, layout = nothing)
    plot_list = []
    for i in 1:n_sim
        model = simple_simulation(environment, sigma, lambda, n_timesteps) 
        push!(plot_list, plot_mixture(model.mixture, environment=environment))
    end
    return plot(plot_list..., size = size, layout=layout)
end

function animate_simulation(simulation_data; layout=nothing, size=nothing, x_limits = [-10, 10],  ymax = .5, show_data = false, speed = 1)
    
    colors = cmap("R1", N=length(simulation_data.update_sequence))
    
    @gif for i in 1:speed:length(simulation_data.update_sequence)
        models = simulation_data.update_sequence[i]
        env = simulation_data.environment
        plot_list = []

        if simulation_data.span[i]
            background = :gray68
        else
            background = :white
        end

        for m in eachindex(models)
            plt = plot_mixture(
                models[m].mixture, 
                title = simulation_data.model_names[m], 
                x_limits=x_limits,
                ymax = ymax, 
                environment = env,
                colors=colors,
                background=background)
            if show_data
                plot!([simulation_data.observations[i]],[0.01],seriestype=:scatter, mc=:black, ms = 5, label=nothing)
                if i > 1 plot!([simulation_data.observations[i-1]],[0.01],seriestype=:scatter, mc=:black, ms = 5, label=nothing, ma=.7) end
                if i > 2 plot!([simulation_data.observations[i-2]],[0.01],seriestype=:scatter, mc=:black, ms = 5, label=nothing, ma=.4) end
                if i > 3 plot!([simulation_data.observations[i-3]],[0.01],seriestype=:scatter, mc=:black, ms = 5, label=nothing, ma=.1) end
            end

            push!(plot_list, plt)
        end
        push!(plot_list,
            plot([2],[2],series_annotations = text("Timestep: "*string(i), :left), legend=false,grid=false,foreground_color_subplot=background)
        )
        push!(plot_list,
            heatmap(hcat(eachindex(simulation_data.update_sequence))', legend=false, color=colors, yaxis=false, xlabel="Color depicting time of creation")
        )
       
        if isnothing(size)
            if isnothing(layout)
                plot(plot_list...)
            else
                plot(plot_list..., layout = layout)
            end
        else
            if isnothing(layout)
                plot(plot_list..., size = size)
            else
                plot(plot_list..., layout = layout, size = size)
            end
        end
    end
    
end

function plot_ambiguity_analysis(df::DataFrame, x_col::Symbol, x_title::String;
                                 metric::Symbol = :ari_mean, line_x = nothing, title = nothing,
                                 clims = startswith(String(metric), "ari") ? (0, 1) : :auto,
                                 reverse_cmap::Bool = false, colorbar::Bool = false,
                                 text_color::Union{Symbol, Nothing} = nothing)
    x_vals = sort(unique(df[!, x_col]))
    dbm_vals = sort(unique(df.DBM))

    plot_data = Matrix{Float64}(undef, length(dbm_vals), length(x_vals))
    for (j, x) in enumerate(x_vals)
        for (i, dbm) in enumerate(dbm_vals)
            row = filter(r -> r[x_col] == x && r.DBM == dbm, df)
            plot_data[i, j] = row[1, metric]
        end
    end

    x_labs = [i % 2 == 0 ? replace(string(x_vals[i]), "-" => "−") : "" for i in eachindex(x_vals)]
    y_labs = string.(dbm_vals)
    maxs_x = [findfirst(x -> x == maximum(row), row) for row in eachrow(plot_data)]

    colorscale = reverse_cmap ? cgrad(cmap("CBL2"), rev=true) : cmap("CBL2")
    main_plot = heatmap(plot_data,
        yticks = (eachindex(y_labs), y_labs), xticks = (eachindex(x_vals), x_labs),
        c = colorscale,
        xlab = x_title,
        ylab = nothing,
        clims = clims,
        left_margin = 1mm,
        bottom_margin = 10mm,
        right_margin = colorbar ? 5mm : -1mm,
        colorbar = colorbar,
        tick_direction = :none,
        xrotation = 30,
        title = title,
        titlefontsize = 10
    )
    if !isnothing(line_x)
        vline!([line_x], color=:white, lw=1, label=nothing, ls=:dash)
    end
    if reverse_cmap
        best_x = [findfirst(x -> x == minimum(row), row) for row in eachrow(plot_data)]
        best_vals = [minimum(plot_data[i,:]) for i in eachindex(best_x)]
    else
        best_x = maxs_x
        best_vals = [maximum(plot_data[i,:]) for i in eachindex(best_x)]
    end
    annotate!([(best_x[i], i,
        Plots.text(
            string(round(best_vals[i], digits=2)),
            8, isnothing(text_color) ? (i > length(dbm_vals) ÷ 2 ? :black : :white) : text_color
        ))
    for i in eachindex(best_x)])

    return main_plot
end

function plot_ambiguity_single(df::DataFrame, x_title::String, metric::Symbol = :ari_mean;
                               clims = (0, 1), reverse_cmap::Bool = false)
    dbm_vals = sort(unique(df.DBM))
    plot_data = [df[df.DBM .== dbm, metric][1] for dbm in dbm_vals]

    colorscale = reverse_cmap ? cgrad(cmap("CBL2"), rev=true) : cmap("CBL2")
    y_labs = string.(dbm_vals)
    side_plot = heatmap(hcat(plot_data),
        yticks = (eachindex(y_labs), y_labs), xticks = ([1], [x_title]),
        c = colorscale,
        xlab = "",
        ylab = "",
        clims = clims,
        left_margin = 0mm,
        bottom_margin = 5mm,
        tick_direction = :none
    )
    annotate!([(1, i,
        Plots.text(
            string(round(plot_data[i], digits=2)),
            8, i > length(dbm_vals) ÷ 2 ? :black : :white
        ))
    for i in eachindex(plot_data)])

    return side_plot
end


function plot_sigma_lambda_analysis(df::DataFrame, metric::Symbol = :ari_mean;
                                    clims = startswith(String(metric), "ari") ? (0, 1) : :auto,
                                    xlab = "λ", ylab = "σ", legend = true, title = "",
                                    reverse_cmap::Bool = false,
                                    mark_baseline::Bool = false,
                                    mark_best::Bool = false)
    sigma_vals = sort(unique(df.sigma))
    lambda_vals = sort(unique(df.lambda))

    plot_data = Matrix{Float64}(undef, length(sigma_vals), length(lambda_vals))
    for (j, λ) in enumerate(lambda_vals)
        for (i, σ) in enumerate(sigma_vals)
            row = filter(r -> r.sigma == σ && r.lambda == λ, df)
            plot_data[i, j] = row[1, metric]
        end
    end

    x_labs = [i % 6 == 2 ? replace(string(lambda_vals[i]), "-" => "−") : "" for i in eachindex(lambda_vals)]
    y_labs = [i % 5 == 2 ? replace(string(sigma_vals[i]), "-" => "−") : "" for i in eachindex(sigma_vals)]

    colorscale = reverse_cmap ? cgrad(cmap("CBL2"), rev=true) : cmap("CBL2")
    plt = heatmap(plot_data,
        yticks = (eachindex(y_labs), y_labs), xticks = (eachindex(lambda_vals), x_labs),
        c = colorscale,
        xlab = xlab,
        ylab = ylab,
        clims = clims,
        legend = legend,
        title = title,
        tick_direction = :none,
        foreground_color_axis = :white,
        foreground_color_border = :white,
    )

    if mark_baseline
        σ_idx = findfirst(v -> isapprox(v, 0.0), sigma_vals)
        λ_idx = findfirst(v -> isapprox(v, 0.7), lambda_vals)
        if !isnothing(σ_idx) && !isnothing(λ_idx)
            scatter!(plt, [λ_idx], [σ_idx], marker=:xcross, color=:red,
                ms=5, msw=2, label=nothing)
        end
    end

    if mark_best
        best_val = reverse_cmap ? minimum(plot_data) : maximum(plot_data)
        idx = findfirst(==(best_val), plot_data)
        scatter!(plt, [idx[2]], [idx[1]], marker=:circle, ms=5,
            markercolor=RGBA(0, 0, 0, 0), markerstrokecolor=:red, markerstrokewidth=2,
            label=nothing)
    end

    return plt
end

function plot_mean_trajectory(df::DataFrame, param::Symbol;
                              metric_mean::Symbol = :ari_mean, baseline_val::Float64 = param == :sigma ? 0.0 : 0.7,
                              ylims = (0, 1), perfect_level = 1, legend = nothing)
    param_vals = sort(unique(df[!, param]), rev=true)
    non_baseline = filter(v -> v != baseline_val, param_vals)
    param_label = param == :sigma ? "σ" : "λ"
    plt = plot(legend = legend)
    lw = 2

    # Plot non-baseline trajectories with colors
    i = 1
    for pval in non_baseline
        sub = sort(filter(r -> r[param] == pval, df), :timestep)
        mean_trajectory = sub[!, metric_mean]
        plot!(mean_trajectory; ylims = ylims, label = replace("$param_label = $pval", "-" => "−"), lw = lw, color = i)
        plot!(mean_trajectory[1:199]; ylims = ylims, label = nothing, lw = lw, color = :gray)
        i += 1
    end

    # Plot baseline as smoothed gray dashed line
    baseline = sort(filter(r -> r[param] == baseline_val, df), :timestep)
    baseline_trajectory = Smoothing.binomial(baseline[!, metric_mean], 10)
    plot!(baseline_trajectory; ylims = ylims, label = nothing, lw = lw, color = :gray, ls = :dash)

    change_start = df.change_start[1]
    change_end = df.change_end[1]
    vspan!([change_start, change_end]; color = :black, label = nothing, alpha = 0.2)
    hline!([perfect_level], color = :black, label = nothing, ls = :dash, lw = 1)
    plt
end

function draw_box!(plt, x, stats; color=:steelblue, width=0.35, lw=1.5)
    hw = width / 2
    # Line from min to max
    plot!(plt, [x, x], [stats.min, stats.max], color=color, lw=lw*0.5, label=nothing)
    # Whisker caps at q01 and q99
    plot!(plt, [x-hw/2, x+hw/2], [stats.q01, stats.q01], color=color, lw=lw, label=nothing)
    plot!(plt, [x-hw/2, x+hw/2], [stats.q99, stats.q99], color=color, lw=lw, label=nothing)
    # Outer box: q05–q95
    plot!(plt, Shape([x-hw, x+hw, x+hw, x-hw], [stats.q05, stats.q05, stats.q95, stats.q95]),
        fillcolor=color, fillalpha=0.05, linecolor=color, lw=lw*0.5, label=nothing)
    # Inner box: q25–q75
    plot!(plt, Shape([x-hw, x+hw, x+hw, x-hw], [stats.q25, stats.q25, stats.q75, stats.q75]),
        fillcolor=color, fillalpha=0.35, linecolor=color, lw=lw, label=nothing)
    # Median line
    plot!(plt, [x-hw, x+hw], [stats.median, stats.median], color=color, lw=lw+1, label=nothing)
    plot!(plt, [x-hw, x+hw], [stats.median, stats.median], color=:black, lw=lw+1, label=nothing, alpha=0.3)
    # Mean diamond
    scatter!(plt, [x-hw], [stats.mean], marker=:diamond, ms=5, color=color, label=nothing)
end

function extract_stats(row, prefix::String)
    (mean   = row[Symbol("$(prefix)_mean")],
     min    = row[Symbol("$(prefix)_min")],
     q01    = row[Symbol("$(prefix)_q01")],
     q05    = row[Symbol("$(prefix)_q05")],
     q25    = row[Symbol("$(prefix)_q25")],
     median = row[Symbol("$(prefix)_median")],
     q75    = row[Symbol("$(prefix)_q75")],
     q95    = row[Symbol("$(prefix)_q95")],
     q99    = row[Symbol("$(prefix)_q99")],
     max    = row[Symbol("$(prefix)_max")])
end

function plot_distribution_comparison(amb_df::DataFrame, baseline_df::Union{DataFrame, Nothing},
                                      perfect_df::Union{DataFrame, Nothing},
                                      metric_prefix::String;
                                      conv_df::Union{DataFrame, Nothing} = nothing,
                                      title="", ylabel="",
                                      show_legend::Bool = true)
    dbm_levels = [6.0, 3.0, 1.0]
    n_groups = length(dbm_levels)

    dbm_labels = ["Unambiguous", "Ambiguous", "Very ambiguous"]

    timesteps = [100, 10000]

    conditions = [
        (label="t = 100",    color=palette(:default)[1]),
        (label="t = 10,000",  color=palette(:default)[2]),
    ]
    if !isnothing(conv_df)
        push!(conditions, (label="t = 100,000", color=palette(:default)[3]))
    end
    if !isnothing(perfect_df)
        push!(conditions, (label="Perfect", color=palette(:default)[4]))
    end
    n_cond = length(conditions)

    group_width = n_cond + 1
    plt = plot(legend=show_legend ? :topright : false, ylabel=ylabel, title=title, grid=:y,
        left_margin=10mm, bottom_margin=5mm, top_margin=5mm, right_margin=7mm,
        tickfontsize=11, guidefontsize=13, titlefontsize=14, legendfontsize=10)

    for (gi, dbm) in enumerate(dbm_levels)
        group_offset = (gi - 1) * group_width
        ci = 0

        row = filter(r -> r.DBM == dbm, amb_df)

        # Boxes for each timestep
        for t in timesteps
            ci += 1
            ts_row = filter(r -> r.n_timesteps == t, row)
            if nrow(ts_row) > 0
                draw_box!(plt, group_offset + ci, extract_stats(ts_row[1, :], metric_prefix), color=conditions[ci].color)
            end

            # Per-box baseline (standard parameterisation) reference line
            if !isnothing(baseline_df)
                base_row = filter(r -> r.DBM == dbm && r.n_timesteps == t, baseline_df)
                if nrow(base_row) > 0
                    base_mean = base_row[1, Symbol("$(metric_prefix)_mean")]
                    x_lo = group_offset + ci - 0.4
                    x_hi = group_offset + ci + 0.4
                    plot!(plt, [x_lo, x_hi], [base_mean, base_mean],
                        color=:black, ls=:dot, lw=3, label=nothing)
                end
            end
        end

        # Convergence (t = 100000)
        if !isnothing(conv_df)
            ci += 1
            conv = filter(r -> r.DBM == dbm, conv_df)
            if nrow(conv) > 0
                draw_box!(plt, group_offset + ci, extract_stats(conv[1, :], metric_prefix), color=conditions[ci].color)
            end
        end

        # Perfect agent
        if !isnothing(perfect_df)
            ci += 1
            perf = filter(r -> r.DBM == dbm, perfect_df)
            if nrow(perf) > 0
                draw_box!(plt, group_offset + ci, extract_stats(perf[1, :], metric_prefix), color=conditions[ci].color)
            end
        end
    end

    # X-axis labels
    group_centers = [(gi - 1) * group_width + (n_cond + 1) / 2 for gi in 1:n_groups]
    plot!(plt, xticks=(group_centers, dbm_labels), xlab="")

    # Legend entries
    for c in conditions
        plot!(plt, [], [], color=c.color, fillalpha=0.4, lw=0, fillrange=0, label=c.label, seriestype=:shape)
    end

    return plt
end

function plot_gmap_comparison(gmap_dfs::Vector, metric_prefix::String;
                              title="", ylabel="", show_legend::Bool = true)
    dbm_levels = [6.0, 3.0, 1.0]
    n_groups = length(dbm_levels)
    dbm_labels = ["Unambiguous", "Ambiguous", "Very ambiguous"]

    n_cond = length(gmap_dfs)
    group_width = n_cond + 1

    plt = plot(legend=show_legend ? :topright : false, ylabel=ylabel, title=title, grid=:y,
        left_margin=10mm, bottom_margin=5mm, top_margin=5mm, right_margin=7mm,
        tickfontsize=11, guidefontsize=13, titlefontsize=14, legendfontsize=10)

    for (gi, dbm) in enumerate(dbm_levels)
        group_offset = (gi - 1) * group_width
        for (ci, gm) in enumerate(gmap_dfs)
            row = filter(r -> r.DBM == dbm, gm.df)
            if nrow(row) > 0
                draw_box!(plt, group_offset + ci, extract_stats(row[1, :], metric_prefix),
                    color=palette(:default)[ci])
            end
        end
    end

    group_centers = [(gi - 1) * group_width + (n_cond + 1) / 2 for gi in 1:n_groups]
    plot!(plt, xticks=(group_centers, dbm_labels), xlab="")

    for (ci, gm) in enumerate(gmap_dfs)
        plot!(plt, [], [], color=palette(:default)[ci], fillalpha=0.4, lw=0,
            fillrange=0, label=gm.label, seriestype=:shape)
    end

    return plt
end

function example_plot(sigma, lambda, n_timesteps, DBM, size; ymax = .5)
    environment = get_simple_environment(Float64(DBM))
    model = init_mixture_model(.01, -6.0, sigma, lambda)
    environment_IDs = eachindex(environment)

    colors = palette(:default)[1:10]

    plots = Vector(undef, n_timesteps+2)
    @inbounds for t in 1:n_timesteps
        i = rand(environment_IDs)
        observation = rand(Normal(environment[i]...))
        update_model!(model, observation)
        p = plot_mixture(model.mixture; environment=environment,ymax = ymax, colors=colors)
        scatter!([observation],[0], c=:black, shape=:vline, msize=10, mswidth=6, bottom_margin = -8mm, top_margin = -2mm)
        hline!([0], c=:black)
        plots[t] = p
    end

    model = simulate_from_model(model, environment, 100-n_timesteps)
    plots[n_timesteps+1] = plot_mixture(model.mixture; environment=environment,ymax = ymax, colors=colors)

    model = simulate_from_model(model, environment, 9900)
    plots[n_timesteps+2] = plot_mixture(model.mixture; environment=environment,ymax = ymax, colors=colors)

    plot(plots..., legend=false, layout=(n_timesteps+2,1), size=size)
end

