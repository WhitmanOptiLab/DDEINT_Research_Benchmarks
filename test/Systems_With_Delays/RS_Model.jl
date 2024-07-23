using DifferentialEquations
using CSV, DataFrames
using Plots
using Dates

ABS_TOL = 1e-9
REL_TOL = 1e-9

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

println("Solving Repressilator model...")
# Solve the DDE problem
# Measured the time taken to solve the problem
start = Dates.now()
sol = solve(prob, MethodOfSteps(Tsit5()), reltol=REL_TOL, abstol=ABS_TOL)
println("It took ", Dates.now() - start, " seconds to solve the Repressilator model")

# Interpolate the solution at finer intervals
t_interp = range(0, stop=50, length=1000)
sol_interp = sol(t_interp)

# Save interpolated solution to CSV
data = DataFrame(time=t_interp, x1=sol_interp[1,:], x2=sol_interp[2,:], x3=sol_interp[3,:])
CSV.write("data/repressilator_model_jl.csv", data)

# print
println("Repressilator model solved successfully")