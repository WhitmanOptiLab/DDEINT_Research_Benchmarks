using DifferentialEquations
using CSV, DataFrames
using Plots

ABS_TOL = 1e-9
REL_TOL = 1e-9

# Breast Cancer Model
function bc_model(du, u, h, p, t)
    p0, q0, v0, d0, p1, q1, v1, d1, d2, beta0, beta1, tau = p
    hist3 = h(p, t - tau)[3]
    du[1] = (v0 / (1 + beta0 * (hist3^2))) * (p0 - q0) * u[1] - d0 * u[1]
    du[2] = (v0 / (1 + beta0 * (hist3^2))) * (1 - p0 + q0) * u[1] +
            (v1 / (1 + beta1 * (hist3^2))) * (p1 - q1) * u[2] - d1 * u[2]
    du[3] = (v1 / (1 + beta1 * (hist3^2))) * (1 - p1 + q1) * u[2] - d2 * u[3]
end

h(p, t) = ones(3)

tau = 1
lags = [tau]

p0 = 0.2;
q0 = 0.3;
v0 = 1;
d0 = 5;
p1 = 0.2;
q1 = 0.3;
v1 = 1;
d1 = 1;
d2 = 1;
beta0 = 1;
beta1 = 1;
p = (p0, q0, v0, d0, p1, q1, v1, d1, d2, beta0, beta1, tau)
tspan = (0.0, 10.0)

u0 = [1.0, 1.0, 1.0]

prob = DDEProblem(bc_model, u0, h, tspan, p; constant_lags = lags)

alg = MethodOfSteps(Tsit5())
# use dormand prince method for stiff 
sol = solve(prob, alg, reltol=REL_TOL, abstol=ABS_TOL)

# Interpolate the solution at finer intervals
t_interp = range(0, stop=10, length=1000)
sol_interp = sol(t_interp)

# Save interpolated solution to CSV
data = DataFrame(time=t_interp, u1=sol_interp[1,:], u2=sol_interp[2,:], u3=sol_interp[3,:])
CSV.write("data/bc_model_jl.csv", data)

# print
println("Breast cancer model solved successfully")

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
sol = solve(prob, MethodOfSteps(Tsit5()), reltol=REL_TOL, abstol=ABS_TOL)

# Interpolate the solution at finer intervals
t_interp = range(0, stop=1000, length=1000)
sol_interp = sol(t_interp)

# Save interpolated solution to CSV
data = DataFrame(time=t_interp, Pa=sol_interp[1,:], Pv=sol_interp[2,:], HR=sol_interp[3,:])
CSV.write("data/cardiovascular_model_jl.csv", data)

# print 
println("Cardiovascular model solved successfully")

# Repressilator Model
tau = 0.1
const beta = 50
const n = 2
const k = 1
const gamma = 1

function repressilator_model(du, u, h, p, t)
    x3_tau = h(p, t - tau)[3]
    x1_tau = h(p, t - tau)[1]
    x2_tau = h(p, t - tau)[2]

    du[1] = beta / (1 + (x3_tau / k)^n) - gamma * u[1]
    du[2] = beta / (1 + (x1_tau / k)^n) - gamma * u[2]
    du[3] = beta / (1 + (x2_tau / k)^n) - gamma * u[3]
end

h(p, t) = ones(3)

u0 = [1.0, 1.0, 1.2]
tspan = (0.0, 50.0)
lags = [tau]

prob = DDEProblem(repressilator_model, u0, h, tspan, nothing; constant_lags=lags)
sol = solve(prob, MethodOfSteps(Tsit5()), reltol=REL_TOL, abstol=ABS_TOL)

# Interpolate the solution at finer intervals
t_interp = range(0, stop=50, length=1000)
sol_interp = sol(t_interp)

# Save interpolated solution to CSV
data = DataFrame(time=t_interp, x1=sol_interp[1,:], x2=sol_interp[2,:], x3=sol_interp[3,:])
CSV.write("data/repressilator_model_jl.csv", data)

# print
println("Repressilator model solved successfully")


# Define the SIR model with delay
# function sir_model(du, u, h, p, t)
#     β, γ, τ = p
#     S_τ = h(p, t - τ)[1]
#     I_τ = h(p, t - τ)[2]
    
#     du[1] = -β * u[1] * u[2]                  # dS/dt
#     du[2] = β * S_τ * I_τ - γ * u[2]           # dI/dt
#     du[3] = γ * u[2]                           # dR/dt
# end

# # Define the history function
# h(p, t) = [1000, 1, 0]  # Initial conditions: S = 1000, I = 1, R = 0

# # Parameters
# β = 0.5
# γ = 0.1
# τ = 1.0
# p = (β, γ, τ)

# # Initial conditions and time span
# u0 = [0.9, 0.1, 0.0]
# tspan = (0.0, 100.0)
# lags = [τ]

# # Define and solve the DDE problem
# prob = DDEProblem(sir_model, u0, h, tspan, p; constant_lags=lags)
# alg = MethodOfSteps(Tsit5())
# sol = solve(prob, alg, reltol=REL_TOL, abstol=ABS_TOL)

# # Interpolate the solution at finer intervals
# t_interp = range(0, stop=100, length=1000)
# sol_interp = sol(t_interp)

# # Save interpolated solution to CSV
# data = DataFrame(time=t_interp, S=sol_interp[1,:], I=sol_interp[2,:], R=sol_interp[3,:])
# CSV.write("data/sir_model_jl.csv", data)

# # print
# println("SIR model with delay solved successfully")