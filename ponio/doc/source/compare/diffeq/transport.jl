using DifferentialEquations
using DelimitedFiles

n_x = 100;
t_span = (0.0, 0.3);
dt = 0.01;

x = 0:(1.0/n_x):1.0;
Δx = x[2] - x[1];

a = 1.0;

function transport!(dy, y, p, t)
  a = p[1]
  Δx = p[2]

  n_x = size(y, 1)

  dy[1] = -a * (y[2] - y[n_x]) / (2.0 * Δx)

  for i = 2:(n_x-1)
    dy[i] = -a * (y[i+1] - y[i-1]) / (2.0 * Δx)
  end

  dy[n_x] = -a * (y[1] - y[n_x-1]) / (2.0 * Δx)
end

y0 = similar(x)
# @. y0 = sin(2.0 * π * x)
for i = 1:n_x
  y0[i] = 0.0
  if 0.25 <= x[i] && x[i] < 0.5
    y0[i] = x[i] - 0.25
  elseif 0.5 <= x[i] && x[i] < 0.75
    y0[i] = -x[i] + 0.75
  end
end

p = [a, Δx]

transport_pb = ODEProblem(transport!, y0, t_span, p)
sol = solve(transport_pb, RK4(), adaptive=false, dt=dt)

y = transpose(mapreduce(transpose, vcat, sol.u))

writedlm("transport.txt", hcat(x, y), " ")
