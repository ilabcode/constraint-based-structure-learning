using Plots; gr()
using Measures
using PerceptualColourMaps
using Statistics
using Smoothing
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

function plot_ambiguity_analysis(results, x_var, x_title, line_x = nothing, title = nothing) # Naming issue with x title and x labels

    plot_data = hcat([results[i] for i in eachindex(x_var)]...)

    x_labs = [i∈2:2:22 ? string(x_var[i]) : "" for i in eachindex(x_var)]
    y_labs = string.(0.5:0.5:6)
    maxs_x = [findfirst(x -> x == maximum(row), row) for row in eachrow(plot_data)]
    main_plot = heatmap(plot_data, 
        yticks = (eachindex(y_labs), y_labs), xticks = (eachindex(x_var), x_labs), 
        c=cmap("CBL2"),
        xlab = x_title,
        ylab = nothing,
        clims = (0, 1),
        left_margin = 1mm,
        bottom_margin = 10mm,
        right_margin = -1mm,
        legend = nothing,
        tick_direction=:none,
        title = title,
        titlefontsize = 10
        )
    if !isnothing(line_x) 
        vline!([line_x], color=:white, lw=1, label=nothing, ls=:dash)
    end
    annotate!([(maxs_x[i], i, 
        Plots.text(
            string(round(maximum(plot_data[i,:]), digits=2))[2:end], 
            8, i>6 ? :black : :white
            )) 
    for i in eachindex(maxs_x)])

    return main_plot
end

function plot_ambiguity_single(results, x_title)

    plot_data = hcat(results)
    y_labs = string.(0.5:0.5:6)
    side_plot = heatmap(plot_data, 
        yticks = (eachindex(y_labs), y_labs), xticks = ([1], [x_title]), 
        c=cmap("CBL2"),
        xlab = "",
        ylab = "",
        clims = (0, 1),
        left_margin = 0mm,
        bottom_margin = 5mm,
        tick_direction=:none
    )
    annotate!([(1, i, 
        Plots.text(
            string(round(results[i], digits=2))[2:end], 
            8, i>5 ? :black : :white
            )) 
        for i in eachindex(results)]) 

    
    return side_plot
end

function plot_sigma_lambda_analysis(results, sigmas, lambdas, max_val = 1; xlab = "λ", ylab = "σ", legend=true, title="")
    plot_data = hcat([results[i] for i in eachindex(lambdas)]...)

    x_labs = [i∈2:6:length(lambdas) ? string(lambdas[i]) : "" for i in eachindex(lambdas)]
    y_labs = [i∈2:5:length(sigmas) ? string(sigmas[i]) : "" for i in eachindex(sigmas)]

    plot = heatmap(plot_data, 
        yticks = (eachindex(y_labs), y_labs), xticks = (eachindex(lambdas), x_labs), 
        c=cmap("CBL2"),
        xlab = xlab,
        ylab = ylab,
        clims = (0, max_val),
        legend=legend,
        title=title,
        )

    return plot
end

function plot_mean_trajectory(results, param, ylims = (0,1), perfect_level = 1, legend=nothing)
    plt = plot(legend=legend)
    lw=2
    i = 1
    for res in results[2:end]
        mat = hcat(res.trajectories...)
        mean_trajectory = vec(mean(mat; dims=2)) # mean across trajectories
        cp = res.change_points[1]
        label = param == "sigma" ? "σ = $(cp.sigma)" : "λ = $(cp.lambda)"
        plot!(mean_trajectory; ylims=ylims, label=label, lw=lw, color = i)
        plot!(mean_trajectory[1:199]; ylims=ylims, label=nothing, lw=lw, color = :gray)
        i += 1
    end

    mat = hcat(results[1].trajectories...)
    mean_trajectory = Smoothing.binomial(vec(mean(mat; dims=2)),10)
    plot!(mean_trajectory; ylims=ylims, label=nothing, lw=lw, color = :gray, ls=:dash)


    cp = results[1].change_points
    vspan!([cp[1].timestep, cp[2].timestep]; color=:black, label=nothing, alpha=.2)
    hline!([perfect_level], color=:black, label=nothing, ls=:dash, lw=1)
    plt
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

