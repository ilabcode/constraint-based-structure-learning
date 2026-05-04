using DataFrames, CSV

function stats_columns!(row::Dict, prefix::String, stats::NamedTuple)
    row["$(prefix)_mean"]   = stats.mean
    row["$(prefix)_min"]    = stats.min
    row["$(prefix)_q01"]    = stats.q01
    row["$(prefix)_q05"]    = stats.q05
    row["$(prefix)_q25"]    = stats.q25
    row["$(prefix)_median"] = stats.median
    row["$(prefix)_q75"]    = stats.q75
    row["$(prefix)_q95"]    = stats.q95
    row["$(prefix)_q99"]    = stats.q99
    row["$(prefix)_max"]    = stats.max
end

function ambiguity_to_dataframe(results::Vector)
    rows = Dict{String,Any}[]
    for res in results
        for (i, dbm) in enumerate(res.DBMs)
            row = Dict{String,Any}(
                "sigma"       => res.sigma,
                "lambda"      => res.lambda,
                "sd"          => res.sd,
                "n_timesteps" => res.n_timesteps,
                "DBM"         => dbm
            )
            stats_columns!(row, "ari",   res.ari[i])
            stats_columns!(row, "mw2",   res.mw2[i])
            stats_columns!(row, "ncomp", res.ncomp[i])
            push!(rows, row)
        end
    end
    return DataFrame(rows)
end

function sigma_lambda_to_dataframe(results::Vector)
    rows = Dict{String,Any}[]
    for res in results
        for (i, sigma) in enumerate(res.sigmas)
            row = Dict{String,Any}(
                "sigma"                 => sigma,
                "lambda"                => res.lambda,
                "distance_between_means"=> res.distance_between_means,
                "n_timesteps"           => res.n_timesteps
            )
            stats_columns!(row, "ari",   res.ari[i])
            stats_columns!(row, "mw2",   res.mw2[i])
            stats_columns!(row, "ncomp", res.ncomp[i])
            push!(rows, row)
        end
    end
    return DataFrame(rows)
end

function trajectory_to_dataframe(results::Vector)
    rows = Dict{String,Any}[]
    for res in results
        cp = res.change_points
        timespan = cp[2].timestep - cp[1].timestep
        T = length(res.ari)
        for t in 1:T
            row = Dict{String,Any}(
                "distance_between_means" => res.distance_between_means,
                "sigma"                  => res.sigma,
                "lambda"                 => res.lambda,
                "timespan"               => timespan,
                "timestep"               => t,
                "change_start"           => cp[1].timestep,
                "change_end"             => cp[2].timestep
            )
            stats_columns!(row, "ari",   res.ari[t])
            stats_columns!(row, "mw2",   res.mw2[t])
            stats_columns!(row, "ncomp", res.ncomp[t])
            push!(rows, row)
        end
    end
    return DataFrame(rows)
end
