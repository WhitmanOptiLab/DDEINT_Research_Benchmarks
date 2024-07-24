module BC_Model

using DifferentialEquations
using CSV, DataFrames
using Plots
using Dates

const ABS_TOL = 1e-9
const REL_TOL = 1e-9

# Breast Cancer Model
@inline function bc_model(du::Vector{Float64}, u::Vector{Float64}, h, p::NTuple{12, Float64}, t::Float64)
    p0, q0, v0, d0, p1, q1, v1, d1, d2, beta0, beta1, tau = p
    hist3 = h(p, t - tau)[3]
    du[1] = (v0 / (1 + beta0 * (hist3^2))) * (p0 - q0) * u[1] - d0 * u[1]
    du[2] = (v0 / (1 + beta0 * (hist3^2))) * (1 - p0 + q0) * u[1] +
            (v1 / (1 + beta1 * (hist3^2))) * (p1 - q1) * u[2] - d1 * u[2]
    du[3] = (v1 / (1 + beta1 * (hist3^2))) * (1 - p1 + q1) * u[2] - d2 * u[3]
end

@inline h(p, t) = ones(3)

const tau = 1.0
const lags = [tau]

const params = (0.2, 0.3, 1.0, 5.0, 0.2, 0.3, 1.0, 1.0, 1.0, 1.0, 1.0, tau)
const tspan = (0.0, 10.0)
const u0 = [1.0, 1.0, 1.0]

function Run_BC_Model()
    prob = DDEProblem(bc_model, u0, h, tspan, params; constant_lags = lags)
    alg = MethodOfSteps(Tsit5())
    
    println("Solving Breast Cancer model...")
    start = Dates.now()
    sol = solve(prob, alg, reltol=REL_TOL, abstol=ABS_TOL)
    elapsed_time = Dates.now() - start
    println("It took ", elapsed_time, " seconds to solve the Breast Cancer model")
    
    t_interp = range(0, stop=10, length=1000)
    sol_interp = sol(t_interp)
    
    data = DataFrame(time=t_interp, u1=sol_interp[1,:], u2=sol_interp[2,:], u3=sol_interp[3,:])
    CSV.write("data/bc_model_jl.csv", data)
    
    println("Breast cancer model solved successfully")
    return elapsed_time
end

end # module
