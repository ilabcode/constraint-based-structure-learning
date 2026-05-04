include("component.jl")
include("mixture.jl")
include("ActionModelSupport.jl")
include("simulation_functions.jl")

using DataFrames, CSV, Plots, StatsPlots

recovery_prior = (
    lambda = Normal(0.7, 2.0),
    sigma  = Normal(0.0, 2.0),
)

lambda_vals = -2.0:0.1:2.0
sigma_vals  = -2.0:0.1:2.0
n_obs_vals  = [100, 400]
dbm_vals    = [3.0, 6.0]

session_key = Symbol("session:single_session")

recovery_results = DataFrame(
    dbm         = Float64[],
    n_obs       = Int[],
    true_lambda = Float64[],
    true_sigma  = Float64[],
    est_lambda  = Float64[],
    est_sigma   = Float64[],
)

for dbm in dbm_vals
    recovery_env = get_simple_environment(dbm)
    generate_observations(n) = [rand(Normal(rand(recovery_env)...)) for _ in 1:n]

    for n_obs in n_obs_vals, λ in lambda_vals, σ in sigma_vals
        sim_model = ActionModel(
            constraint_based_structure_learning,
            submodel = MixModelSubmodel(sigma = σ, lambda = λ),
            parameters = (; sigma = Parameter(σ), lambda = Parameter(λ)),
            observations = (; observation = Observation()),
            actions = (; report = Action(Categorical)),
        )
        sim_agent = init_agent(sim_model)
        obs     = generate_observations(n_obs)
        actions = simulate!(sim_agent, obs)

        modelfit = create_model(sim_model, recovery_prior, obs, actions)
        sample_posterior!(modelfit; n_chains = 1, n_samples = 2000, sampler = MH(), init_params = nothing)

        sp = get_session_parameters!(modelfit)
        est_λ = median(getproperty(sp.value.lambda, session_key))
        est_σ = median(getproperty(sp.value.sigma,  session_key))

        push!(recovery_results, (dbm, n_obs, λ, σ, est_λ, est_σ))
    end
end

# CSV.write("data/recovery_results.csv", recovery_results)
recovery_results = CSV.read("data/recovery_results.csv", DataFrame)

function recovery_boxplot(df, param::Symbol, dbm::Float64, n_obs::Int, vals)
    sub = filter([:dbm, :n_obs] => (d, n) -> d == dbm && n == n_obs, df)
    true_col = Symbol("true_", param)
    est_col  = Symbol("est_",  param)
    sym      = param == :lambda ? "λ" : "σ"
    p = @df sub boxplot(cols(true_col), cols(est_col),
        label = nothing,
        xlabel = "True $sym", ylabel = "Estimated $sym",
        title = "$(uppercase(string(param))) - DBM=$dbm, N=$n_obs",
        bar_width = 0.1, c=:honeydew2)
    plot!(vals, collect(vals), c = :maroon, linewidth = 4, label = nothing, alpha = 0.5)
    p
end

plots = [
    recovery_boxplot(recovery_results, :lambda, 3.0, 100, lambda_vals)  recovery_boxplot(recovery_results, :sigma, 3.0, 100, sigma_vals);
    recovery_boxplot(recovery_results, :lambda, 3.0, 400, lambda_vals)  recovery_boxplot(recovery_results, :sigma, 3.0, 400, sigma_vals);
    recovery_boxplot(recovery_results, :lambda, 6.0, 100, lambda_vals)  recovery_boxplot(recovery_results, :sigma, 6.0, 100, sigma_vals);
    recovery_boxplot(recovery_results, :lambda, 6.0, 400, lambda_vals)  recovery_boxplot(recovery_results, :sigma, 6.0, 400, sigma_vals)
]

plot(plots..., layout = (4, 2), size = (1200, 1800))
