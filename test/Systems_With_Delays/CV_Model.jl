module CV_Model

using DifferentialEquations
using CSV, DataFrames
using Plots
using Dates

const ABS_TOL = 1e-9
const REL_TOL = 1e-9

# Cardiovascular Model
struct Params
    ca::Float64
    cv::Float64
    R::Float64
    r::Float64
    Vstr::Float64
    alpha0::Float64
    alphas::Float64
    alphap::Float64
    alphaH::Float64
    beta0::Float64
    betas::Float64
    betap::Float64
    betaH::Float64
    gammaH::Float64
end

const params = Params(1.55, 519, 1.05, 0.068, 67.9, 93, 93, 93, 0.84, 7, 7, 7, 1.17, 0)
const tau = 4.0

@inline function ddefun(du::Vector{Float64}, u::Vector{Float64}, h, p::Params, t::Float64)
    R = t <= 600 ? 1.05 : 0.21 * exp(600 - t) + 0.84
    ylag = h(p, t - tau)
    Patau = ylag[1]
    Paoft = u[1]
    Pvoft = u[2]
    Hoft = u[3]

    du[1] = - (1 / (p.ca * R)) * Paoft + (1 / (p.ca * R)) * Pvoft + (1 / p.ca) * p.Vstr * Hoft
    du[2] = (1 / (p.cv * R)) * Paoft - (1 / (p.cv * R) + 1 / (p.cv * p.r)) * Pvoft
    Ts = 1 / (1 + (Patau / p.alphas)^p.betas)
    Tp = 1 / (1 + (p.alphap / Paoft)^p.betap)
    du[3] = (p.alphaH * Ts) / (1 + p.gammaH * Tp) - p.betaH * Tp
end

@inline function history(p::Params, t::Float64)
    P0 = 93
    Paval = P0
    Pvval = (1 / (1 + p.R / p.r)) * P0
    Hval = (1 / (p.R * p.Vstr)) * (1 / (1 + p.r / p.R)) * P0
    return [Paval, Pvval, Hval]
end

const u0 = [93, (1 / (1 + params.R / params.r)) * 93, (1 / (params.R * params.Vstr)) * (1 / (1 + params.r / params.R)) * 93]
const tspan = (0.0, 1000.0)

function Run_CV_Model()
    prob = DDEProblem(ddefun, u0, history, tspan, params; constant_lags = [tau], tstops = [600])
    alg = MethodOfSteps(Tsit5())

    println("Solving Cardiovascular model...")
    start = Dates.now()
    sol = solve(prob, alg, reltol=REL_TOL, abstol=ABS_TOL)
    elapsed_time = Dates.now() - start
    println("It took ", elapsed_time, " seconds to solve the Cardiovascular model")
    
    t_interp = range(0, stop=1000, length=1000)
    sol_interp = sol(t_interp)
    
    data = DataFrame(time=t_interp, Pa=sol_interp[1,:], Pv=sol_interp[2,:], HR=sol_interp[3,:])
    CSV.write("data/cardiovascular_model_jl.csv", data)
    
    println("Cardiovascular model solved successfully")
    return elapsed_time
end

end # module
