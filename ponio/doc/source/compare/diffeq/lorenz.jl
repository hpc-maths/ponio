using DifferentialEquations
using DelimitedFiles

function lorenz!(dy, y, p, t)
  σ = p[1]
  ρ = p[2]
  β = p[3]

  dy[1] = σ * (y[2] - y[1])
  dy[2] = y[1] * (ρ - y[3]) - y[2]
  dy[3] = y[1] * y[2] - β * y[3]
end

y0 = [1.0; 1.0; 1.0]
p = [10.0, 28.0, 8.0 / 3.0]
t_span = (0.0, 20.0)
dt = 0.05

lorenz_pb = ODEProblem(lorenz!, y0, t_span, p)
sol = solve(lorenz_pb, RK4(), adaptive=false, dt=dt)

y = mapreduce(permutedims, vcat, sol.u)

writedlm("lorenz.txt", hcat(sol.t, y), " ")
