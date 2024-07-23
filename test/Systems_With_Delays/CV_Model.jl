using DifferentialEquations
using CSV, DataFrames
using Plots
using Dates

ABS_TOL = 1e-9
REL_TOL = 1e-9

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

p = Params(1.55, 519, 1.05, 0.068, 67.9, 93, 93, 93, 0.84, 7, 7, 7, 1.17, 0)

tau = 4.0

function ddefun(du, u, h, p, t)
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

function history(p, t)
    P0 = 93
    Paval = P0
    Pvval = (1 / (1 + p.R / p.r)) * P0
    Hval = (1 / (p.R * p.Vstr)) * (1 / (1 + p.r / p.R)) * P0
    return [Paval, Pvval, Hval]
end

u0 = [93, (1 / (1 + p.R / p.r)) * 93, (1 / (p.R * p.Vstr)) * (1 / (1 + p.r / p.R)) * 93]
tspan = (0.0, 1000.0)

prob = DDEProblem(ddefun, u0, history, tspan, p; constant_lags = [tau], tstops = [600])

println("Solving Cardiovascular model...")

start = Dates.now()
sol = solve(prob, MethodOfSteps(Tsit5()), reltol=REL_TOL, abstol=ABS_TOL)
println("It took ", Dates.now() - start, " seconds to solve the Cardiovascular model")

# Interpolate the solution at finer intervals
t_interp = range(0, stop=1000, length=1000)
sol_interp = sol(t_interp)

# Save interpolated solution to CSV
data = DataFrame(time=t_interp, Pa=sol_interp[1,:], Pv=sol_interp[2,:], HR=sol_interp[3,:])
CSV.write("data/cardiovascular_model_jl.csv", data)

# print 
println("Cardiovascular model solved successfully")
