using DifferentialEquations ,  StochasticDiffEq, DiffEqJump, Plots
###------------------
g = 0.05
function INPLACE_WHITE_NOISE_DIST_Trunc(rand_vec,W,dt,rng)
    rand_vec_1 = (g * dt).* rand_vec
  DiffEqNoiseProcess.wiener_randn!(rng,rand_vec)
  rand_vec .*= sqrt(abs(dt))
    rand_vec .= min.(rand_vec_1, rand_vec)
end
TruncatedWienerProcess(t0,W0,Z0=nothing) = NoiseProcess(t0,W0,Z0,INPLACE_WHITE_NOISE_DIST_Trunc,
    nothing)
###------------------
struct AffectIndex{F1, F2}
affect_index!::F1
index::F2
end
function (p::AffectIndex)(integrator)
return p.affect_index!(integrator, p.index)
end

####
function μ_SDE(du,u,p,t)
du .= p.μ
end

function σ_SDE(du,u,p,t)
du .= p.σ
end
p = (μ = 0.01, σ = 0.1, N = 5) # if all constant
T = 10.0  # maximum time length
x_iv = rand(p.N)  # just draws from the inital condition
#with TruncatedWienerProcess
w = TruncatedWienerProcess(0.0,zeros(p.N),zeros(p.N))
prob = SDEProblem(μ_SDE, σ_SDE, x_iv ,(0.0, T), p, noise = w)
#without TruncatedWienerProcess
#prob = SDEProblem(μ_SDE, σ_SDE, x_iv ,(0.0, T), p)


rate(u,p,t) = 0.2
affect_index!(integrator, index) = (integrator.u[index] =
    max(integrator.u[index], integrator.u[rand(1:integrator.p.N)]) )
jumps = [ConstantRateJump(rate,AffectIndex(affect_index!, i)) for i in 1:p.N]
jump_prob = JumpProblem(prob,DirectFW(),JumpSet((), jumps, nothing, nothing))
sol = solve(jump_prob, SRIW1(), adaptive=false, dt = 0.1);
#sol = solve(jump_prob, SRIW1());#without TruncatedWienerProcess
plot(sol)
