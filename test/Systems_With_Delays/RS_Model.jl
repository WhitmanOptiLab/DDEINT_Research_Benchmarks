module RS_Model

using DifferentialEquations
using CSV, DataFrames
using Dates
using Plots

const ABS_TOL = 1e-9
const REL_TOL = 1e-9

# Repressilator Model
const tau = 0.1
const beta = 50.0
const n = 2.0
const k = 1.0
const gamma = 1.0

@inline function repressilator_model(du::Vector{Float64}, u::Vector{Float64}, h, p, t::Float64)
    x3_tau = h(p, t - tau)[3]
    x1_tau = h(p, t - tau)[1]
    x2_tau = h(p, t - tau)[2]

    du[1] = beta / (1 + (x3_tau / k)^n) - gamma * u[1]
    du[2] = beta / (1 + (x1_tau / k)^n) - gamma * u[2]
    du[3] = beta / (1 + (x2_tau / k)^n) - gamma * u[3]
end

@inline h(p, t) = ones(3)

const u0 = [1.0, 1.0, 1.2]
const tspan = (0.0, 50.0)
const lags = [tau]

function Run_Repressilator_Model()
    prob = DDEProblem(repressilator_model, u0, h, tspan, nothing; constant_lags=lags)

    println("Solving Repressilator model...")
    start = Dates.now()
    sol = solve(prob, MethodOfSteps(Tsit5()), reltol=REL_TOL, abstol=ABS_TOL)
    elapsed_time = Dates.now() - start
    println("It took ", elapsed_time, " seconds to solve the Repressilator model")
    
    t_interp = range(0, stop=50, length=1000)
    sol_interp = sol(t_interp)
    
    data = DataFrame(time=t_interp, x1=sol_interp[1,:], x2=sol_interp[2,:], x3=sol_interp[3,:])
    CSV.write("data/repressilator_model_jl.csv", data)
    
    println("Repressilator model solved successfully")
    return elapsed_time
end

end # module
