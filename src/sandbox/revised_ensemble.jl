using DifferentialEquations, Plots
using DifferentialEquations.EnsembleAnalysis  # Needed? just to make sure we have acces...
function μ_SDE(du,u,p,t)
  du .= p.μ
end

function σ_SDE(du,u,p,t)
  du .= p.σ
end

p = (μ = 0.01, σ = 0.1, N = 3) 
T = 10.0
x_iv = rand(p.N)dd

prob = SDEProblem(μ_SDE, σ_SDE, x_iv ,(0.0, T), p)
function output_func(sol, i)
    sol.u .= [row[1:2] for row in sol.u]
    return (sol, false)
end

ensemble_prob = EnsembleProblem(prob, output_func = output_func)
sim = solve(ensemble_prob,Tsit5(), EnsembleSerial(), trajectories = 2)

#Extract results
summ = EnsembleSummary(sim, quantiles=[0.05, 0.95])
summ.u.u # mean?
