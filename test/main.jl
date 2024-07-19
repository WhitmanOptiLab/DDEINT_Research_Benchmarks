using DifferentialEquations
using CSV, DataFrames
using Plots

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
sol = solve(prob, alg)

# save sol to a csv file
CSV.write("data/bc_model_jl.csv", DataFrame(sol))

# plot the solution
# plot(sol, vars=(0, 1), label="u1")
# plot!(sol, vars=(0, 2), label="u2")
# plot!(sol, vars=(0, 3), label="u3")

# # save the plot
# savefig("plots/bc_model.png")

# Define physical parameters
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

# Define delay
tau = 4.0

# Define the DDE function
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

# Define the history function
function history(p, t)
    P0 = 93
    Paval = P0
    Pvval = (1 / (1 + p.R / p.r)) * P0
    Hval = (1 / (p.R * p.Vstr)) * (1 / (1 + p.r / p.R)) * P0
    return [Paval, Pvval, Hval]
end

# Initial conditions and time span
u0 = [93, (1 / (1 + p.R / p.r)) * 93, (1 / (p.R * p.Vstr)) * (1 / (1 + p.r / p.R)) * 93]
tspan = (0.0, 1000.0)

# Define and solve the DDE problem
prob = DDEProblem(ddefun, u0, history, tspan, p; constant_lags = [tau], tstops = [600])
sol = solve(prob, MethodOfSteps(Tsit5()))

# save sol to a csv file
CSV.write("data/cardiovascular_model_jl.csv", DataFrame(sol))

# # Plot the solution
# plot(sol, vars=(0, 3), label="H(t)", xlabel="Time t", ylabel="Heart Rate")
# title!("Heart Rate for Baroreflex-Feedback Mechanism")

# # Save the plot to a file
# savefig("plots/cardiovascular_model.png")
