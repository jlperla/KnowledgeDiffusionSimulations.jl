using DifferentialEquations ,BenchmarkTools,  StochasticDiffEq, DiffEqJump, Plots
struct AffectIndex{F1, F2}
  affect_index!::F1
  index::F2
end
function (p::AffectIndex)(integrator)
    return p.affect_index!(integrator, p.index)
end

####
function affect!(integrator)
  N = length(integrator.u)
  n = rand(1:N)
  n2 = rand(1:N)
  integrator.u[n] = max(integrator.u[n], integrator.u[n2])
end

function μ_SDE(du,u,p,t)
  du .= p.μ
end

function σ_SDE(du,u,p,t)
  du .= p.σ
end

function test_many_jumps()
    p = (μ = 0.01, σ = 0.1, N = 10) # if all constant
    T = 10.0  # maximum time length
    x_iv = rand(p.N)  # just draws from the inital condition

    prob = SDEProblem(μ_SDE, σ_SDE, x_iv ,(0.0, T), p)
    rate(u,p,t) = 0.2
    affect_index!(integrator, index) = (integrator.u[index] =
                                    max(integrator.u[index], integrator.u[rand(1:integrator.p.N)]) )
    jumps = [ConstantRateJump(rate,AffectIndex(affect_index!, i)) for i in 1:p.N]
    jump_prob = JumpProblem(prob,DirectFW(),JumpSet((), jumps, nothing, nothing))
    @btime solve($jump_prob, SRIW1());
end

function test_single_jump()
    p = (μ = 0.01, σ = 0.1, N = 10) # if all constant
    T = 10.0  # maximum time length
    x_iv = rand(p.N)  # just draws from the inital condition

    prob = SDEProblem(μ_SDE, σ_SDE, x_iv ,(0.0, T), p)
    rate(u,p,t) = 0.2*p.N
    jump = ConstantRateJump(rate,affect!)
    jump_prob = JumpProblem(prob,DirectFW(),jump)
    @btime solve($jump_prob, SRIW1());
end
test_many_jumps();
test_single_jump();
