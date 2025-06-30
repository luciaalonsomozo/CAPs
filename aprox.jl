# aproximacion.jl

using Plots
using RadiiPolynomial
using LinearAlgebra


f(x, a) = x^2 + a# ecuación diferencial 1
Df(x) = ExactReal(2)*x
∂af = [1]

# u = (x,a) para la f
F(u, v, w) = [f(u[1],u[2]); (u - w)⋅v]
DF(u, v) = [Df(u[1]) ∂af; transpose(v)]

# initial point on the branch of equilibria

a_init = 0.2

u_init = 0
u_init, success = newton(u -> (f(u, a_init), Df(u)), u_init)

# next point on the branch of equilibria

x_init = [u_init]
v = vec(nullspace([Df(x_init[1]) ∂af]))[:, 1]  # Take the first column if it's a matrix


println("Size of x_init: ", size(x_init))
println("Size of v: ", size(v))


δ = 5e-2 # step size

w = x_init + δ * v # predictor

x_final, success = newton(x -> (F(x, v, w), DF(x, v)), w)

# Arc-length parametrization: setting up chebyshev nodes for smooth interpolation of curve
N = 700
N_fft = nextpow(2, 2N + 1)
npts = N_fft ÷ 2 + 1
arclength = 15.0
arclength_grid = [0.5 * arclength - 0.5 * cospi(2j/N_fft) * arclength for j ∈ 0:npts-1]

# x_grid = List of points along the solution curve, v_grid = Tangent vectors (directions of continuation)
x_grid = Vector{Vector{Float64}}(undef, npts)
v_grid = Vector{Vector{Float64}}(undef, npts)

# initialize: first point x_init and initial tangent computed via nullspace

direction = [0, 0, -1] # starts by decreasing the parameter
x_grid[1] = x_init
v_grid[1] = vec(nullspace([Df(x_grid[1][1]) ∂af]))
if direction ⋅ v_grid[1] < 0 # enforce direction
    v_grid[1] .*= -1
end

# run continuation scheme: The predictor step is simple Euler step, and the correction step is predicted via Newton

for i ∈ 2:npts
    δᵢ = arclength_grid[i] - arclength_grid[i-1]

    wᵢ = x_grid[i-1] .+ δᵢ .* v_grid[i-1]

    x, success = newton(x -> (F(x, v_grid[i-1], wᵢ), DF(x, v_grid[i-1])), wᵢ; verbose = true)
    success || error()

    x_grid[i] = x
    v_grid[i] = vec(nullspace([Df(x_grid[i][1]) ∂af]))
    if v_grid[i-1] ⋅ v_grid[i] < 0 # keep the same direction
        v_grid[i] .*= -1
    end
end

# construct the approximations

grid2cheb(x_fft::Vector{<:Vector}, N) =
    [rifft!(complex.(getindex.(x_fft, i)), Chebyshev(N)) for i ∈ eachindex(x_fft[1])]

grid2cheb(x_fft::Vector{<:Matrix}, N) =
    [rifft!(complex.(getindex.(x_fft, i, j)), Chebyshev(N)) for i ∈ axes(x_fft[1], 1), j ∈ axes(x_fft[1], 2)]

# Chebyshev interpolation. Extends grid for symmetry. 
x_fft = [reverse(x_grid) ; x_grid[begin+1:end-1]]
x̄ = map(x -> interval.(x), grid2cheb(x_fft, N))

v_fft = [reverse(v_grid) ; v_grid[begin+1:end-1]]
v̄ = map(v -> interval.(v), grid2cheb(v_fft, N))

# Approximation of inverse Jacobians along the curve
A = map(A -> interval.(A), grid2cheb(inv.(DF.(x_fft, v_fft)), N))

# compute the bounds

function cheb2grid(x::VecOrMat{<:Sequence}, N_fft)
    vals = fft.(x, N_fft)
    return [real.(getindex.(vals, i)) for i ∈ eachindex(vals[1])]
end

# AF is a polynomial with respect to s of order 4N

N4 = 4N
N4_fft = nextpow(2, 2N4 + 1)

AF_fft = cheb2grid(A, N4_fft) .* F.(cheb2grid(x̄, N4_fft), cheb2grid(v̄, N4_fft), cheb2grid(x̄, N4_fft))
AF = grid2cheb(AF_fft, N4)

Y = norm(norm.(AF, 1), 1)

# ADF is a polynomial with respect to s of order 3N

N3 = 3N
N3_fft = nextpow(2, 2N3 + 1)

I_ADF_fft = [I] .- cheb2grid(A, N3_fft) .* DF.(cheb2grid(x̄, N3_fft), cheb2grid(v̄, N3_fft))
I_ADF = grid2cheb(I_ADF_fft, N3)

Z₁ = opnorm(norm.(I_ADF, 1), 1)

R = 1.2sup(Y)

a, ϵ = 5, 1
Z₂ = opnorm(norm.(A, 1), 1) * max(abs(2a + 2) + 6(norm(x̄[1], 1) + R) + abs(ϵ), abs(ϵ))

setdisplay(:full)

interval_of_existence(Y, Z₁, Z₂, R)