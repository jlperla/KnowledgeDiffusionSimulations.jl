using DifferentialEquations, StochasticDiffEq, DiffEqJump, Plots, DiffEqNoiseProcess
function INPLACE_WHITE_NOISE_DIST_Trunc(g, rand_vec,W,dt,rng)
  DiffEqNoiseProcess.wiener_randn!(rng,rand_vec)
  rand_vec .*= sqrt(abs(dt))
    rand_vec .= min.(rand_vec, g * abs(dt))
end
TruncatedWienerProcess(g, t0,W0, Z0=nothing) = NoiseProcess(t0,W0,Z0,
    (rand_vec,W,dt,rng) -> INPLACE_WHITE_NOISE_DIST_Trunc(g, rand_vec,W,dt,rng),nothing)

function μ_SDE(u,p,t)
    p.μ
end

function σ_SDE(u,p,t)
    p.σ
end  # just draws from the inital condition

p = (μ = 0.01, σ = 0.1, N = 5, g = 0.05) # if all constant
T = 10.0
x_iv = rand(p.N)
w = TruncatedWienerProcess(p.g, 0.0, zeros(p.N), zeros(p.N))
prob = SDEProblem(μ_SDE, σ_SDE, x_iv ,(0.0, T), p, noise = w)
sol = solve(prob, EM(), dt= 0.1)
plot(sol, label = "g = 0.05")
