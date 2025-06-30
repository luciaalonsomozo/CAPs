using Plots
using RadiiPolynomial

using PlutoTeachingTools, PlutoUI # packages for the notebook

function Λ⁻¹(β, space)
	∂ = project(Derivative(1), space, space) 
	∂² = ∂ * ∂
	return -inv(∂² + β[1] * ∂ + β[2]*I)
end

F(x, β, c, space) = x + Λ⁻¹(β, space) * (x*x + β[3] * c)

DF(x, β, domain, codomain) = I + Λ⁻¹(β, codomain) * project(Multiplication(2x), domain, codomain)

β = [0.1, 4.0, 1.0]
c = Sequence(Fourier(1, 1.0), [0.5, 0.0, 0.5])
K = 10
initial_guess = zeros(ComplexF64, Fourier(K, 1.0));
F_DF(x) = (F(x, β, c, space(x)), DF(x, β, space(x), space(x)))
x_bar, success = newton(F_DF, initial_guess)

println("x_bar = ", x_bar)

begin
	plot(LinRange(0, 2π, 201), t -> real(x_bar(t)); linewidth = 3, label = "x_bar")
	xlims!(0, 2π)
	xlabel!("t")
end

A_K = inv(DF(x_bar, β, space(x_bar), space(x_bar)))
conjugacy_symmetry!(x_bar)
conj(x_bar[K:-1:-K]) == x_bar[-K:K]
iβ = [interval(1)/interval(10), interval(4), interval(1)]

ic = interval(c)
ix_bar = interval(x_bar)
iA_K = interval(A_K)
ν = interval(1.1) # does not have to be exactly 11/10
X = Ell1(GeometricWeight(ν))

function bound_Y(x_bar, β, c, A, X)
	K = order(x_bar)
	space = Fourier(2K, frequency(x_bar))

	Fx = F(x_bar, β, c, space)

	Id = interval(I)
	ext_A = complex(project(Id, space, space))
	ext_A[-K:K,-K:K] = coefficients(A)

	return norm(ext_A * Fx, X)
end

Y = bound_Y(ix_bar, iβ, ic, iA_K, X)

function bound_Z₁(x_bar, β, A, X)
	K = order(x_bar)
	domain = Fourier(2K, frequency(x_bar))
	codomain = Fourier(3K, frequency(x_bar))

	DFx = DF(x_bar, β, domain, codomain)

	Id = interval(I)
	ext_A = complex(project(Id, codomain, codomain))
	ext_A[-K:K,-K:K] = coefficients(A)

	return opnorm(Id - ext_A * DFx, X)
end

Z₁ = bound_Z₁(ix_bar, iβ, iA_K, X)

R = Inf
function bound_Z₂(β, A, X)
	K = order(domain(A))
	return 2max(opnorm(Λ⁻¹(β, domain(A)) * A, X), abs((K+1)^2 - im * β[1] * (K+1) - β[2]))
end

Z₂ = bound_Z₂(iβ, iA_K, X)

ie = interval_of_existence(Y, Z₁, Z₂, R)
r = inf(ie)
println("r = ", r)

rigorous_eval(x, r, t) = interval(real(x(t)), r; format = :midpoint)

ts = LinRange(0, 2π, 101)

display(plot([[interval(ts[i], ts[i+1]), rigorous_eval(ix_bar, r, interval(ts[i], ts[i+1]))] for i = 1:length(ts)-1];
	label = "", color = :palegreen3))

plot!(ts, t -> real(x_bar(t));
	label = "x_bar", color = :royalblue1, linewidth = 3)

xlims!(0, 2π)
xlabel!("t")

sleep(10)
