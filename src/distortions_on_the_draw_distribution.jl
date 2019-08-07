using DifferentialEquations ,BenchmarkTools,  StochasticDiffEq, DiffEqJump, Plots
struct AffectIndex{F1, F2}
  affect_index!::F1
  index::F2
end
function (p::AffectIndex)(integrator)
    return p.affect_index!(integrator, p.index)
end

####
struct RateIndex{F1, F2}
  rate_index::F1
  index::F2
end
function (q::RateIndex)(u, p, t)
    return q.rate_index(u, p, t, q.index)
end
function μ_SDE(du,u,p,t)
  du .= p.μ
end

function σ_SDE(du,u,p,t)
  du .= p.σ
end

function affect!(integrator)
    κ = 0.5
    N = length(integrator.u)
    sort_u = sort(integrator.u)
    F = range(0.0, 1, length = N)
    draw_quantile = rand(2)  # draw quantile on (0,1)
    u_n = sort_u[findlast(q -> q <= draw_quantile[1], F.^κ)]
    u_n_2 = sort_u[findlast(q -> q <= draw_quantile[2], F.^κ)]
    integrator.u[integrator.u .== u_n] =
                max(integrator.u[integrator.u .== u_n], integrator.u[integrator.u .== u_n_2])
end

p = (μ = 0.01, σ = 0.1, N = 5) # if all constant
T = 10.0  # maximum time length
x_iv = rand(p.N)  # just draws from the inital condition

prob = SDEProblem(μ_SDE, σ_SDE, x_iv ,(0.0, T), p)
rate(u,p,t) = 0.2*p.N
jump = ConstantRateJump(rate,affect!)
jump_prob = JumpProblem(prob,DirectFW(),jump)
sol = solve(jump_prob, SRIW1());
plot(sol)
