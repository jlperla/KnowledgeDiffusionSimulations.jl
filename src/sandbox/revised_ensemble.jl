using DifferentialEquations, Plots
using DifferentialEquations.EnsembleAnalysis
using StatsBase, DataFrames

# DE setup
function μ_SDE(du,u,p,t)
  du .= p.μ
end

function σ_SDE(du,u,p,t)
  du .= p.σ
end

# Callback setup
T = 10.;
function save_func(u, t, integrator) 
    if length(integrator.p.moments) == 0 
        g = 0.
    else
        g = mean(u) - integrator.p.moments[end][2]
    end            
    moments = [minimum(u), mean(u), maximum(u), g]
    push!(integrator.p.moments, moments) 
end 

p = (μ = 0.01, σ = 0.1, N = 3, moments = Array{Array{Float64, 1}, 1}());
x_iv = rand(p.N)
prob = SDEProblem(μ_SDE, σ_SDE, x_iv ,(0.0, T), p);
saveat = 0:0.1:T

cb = FunctionCallingCallback(save_func;
                 funcat=saveat,
                 func_everystep=true,
                 func_start = true,
                 tdir=1);
    
# Ensemble setup 
function output_func(sol, i)
    return (u = sol.prob.p.moments, t = sol.t), false
end
    
function prob_func(prob,i,repeat)
    p = (μ = 0.01, σ = 0.1, N = 3, moments = Array{Array{Float64, 1}, 1}());
    SDEProblem(μ_SDE, σ_SDE, x_iv ,(0.0, T), p);
end
    
ensemble_prob = EnsembleProblem(prob, prob_func = prob_func, output_func = output_func)
sim = solve(ensemble_prob,Tsit5(), EnsembleSerial(), trajectories = 2, callback = cb)