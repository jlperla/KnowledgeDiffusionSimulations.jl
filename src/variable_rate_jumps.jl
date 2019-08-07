using DifferentialEquations ,BenchmarkTools,  StochasticDiffEq, DiffEqJump, Plots
struct AffectIndex{F1, F2}
  affect_index!::F1
  index::F2
end
function (p::AffectIndex)(integrator)
    return p.affect_index!(integrator, p.index)
end
###
struct RateIndex{F1, F2}
  rate_index::F1
  index::F2
end
function (q::RateIndex)(u, p, t)
    return q.rate_index(u, p, t, q.index)
end
###
function rate_index(u, p, t, index)
   u_max = maximum(u)
   u_min = minimum(u)
   return p.ρ_max* ((u[index] - u_max)^2 / (u_max - u_min))
end
###
affect_index!(integrator, index) = (integrator.u[index] =
    max(integrator.u[index], integrator.u[rand(1:integrator.p.N)]) )
###
function μ_SDE(du,u,p,t)
  du .= p.μ
end

function σ_SDE(du,u,p,t)
  du .= p.σ
end
p = (μ = 0.01, σ = 0.1, N = 2, ρ_max = 2.0) # if all constant
T = 10.0  # maximum time length
x_iv = rand(p.N)  # just draws from the inital condition
prob = SDEProblem(μ_SDE, σ_SDE, x_iv ,(0.0, T), p)

jumps = [VariableRateJump(RateIndex(rate_index, i), AffectIndex(affect_index!, i)) for i in 1:p.N];
jump_prob = JumpProblem(prob,DirectFW(),JumpSet((jumps),(),nothing,nothing))
sol = solve(jump_prob, SRIW1());
sol_hcat = hcat(sol.u...)
plot(sol.t, transpose(sol_hcat[1:2, :]), label = "Single Jump", legend = false)
