# orbita3.jl

using Plots
using RadiiPolynomial

include("estabilidad.jl")  # Incluir el archivo
using .Estabilidad  # Importar el módulo

f(x, μ) = μ * x * (1 - x) # logistic map

function F(x, μ)
    v = zeros(eltype(x), 3)
    v[1] = f(x[1], μ) - x[2]
	v[2] = f(x[2], μ) - x[3]
	v[3] = f(x[3], μ) - x[1]
	return v
end

Df(x, μ) = μ * (1 - 2x) # derivative of the logistic map w.r.t. x

function DF(x, μ) # derivative of the zero-finding problem w.r.t. x
	M = zeros(eltype(x), 3, 3)
	M[1,1] = Df(x[1], μ)
	M[1,2] = -1
	M[2,2] = Df(x[2], μ)
	M[2,3] = -1
	M[3,1] = -1
	M[3,3] = Df(x[3], μ)
	return M
end

initial_guess = [0.5, 1.0, 0.1]
μ  = 3.9
F_DF(x) = (F(x, μ), DF(x, μ)) # returns `F` and `DF` at `μ = 3.9`
x_bar, success = newton(F_DF, initial_guess)
print(x_bar)

A = inv(DF(x_bar, μ))
# since 39/10 is not representable as a floating-point numbercomo 
iμ = interval(39)/interval(10)
ix_bar = interval(x_bar)
iA = interval(A)
Y = norm(iA * F(ix_bar, iμ), 1)
Z₁ = opnorm(I - iA * DF(ix_bar, iμ), 1)
R = Inf
Z₂ = interval(2) * iμ * opnorm(iA, 1)
ie = interval_of_existence(Y, Z₁, Z₂, R)
r = inf(ie)

println("r = ", r)

enclosure_x_tilde_1 = interval(ix_bar[1], r; format = :midpoint)
enclosure_x_tilde_2 = interval(ix_bar[2], r; format = :midpoint)
enclosure_x_tilde_3 = interval(ix_bar[3], r; format = :midpoint)

if isdisjoint_interval(enclosure_x_tilde_1, enclosure_x_tilde_2, enclosure_x_tilde_3)
	println("Orbita de periodo 3 encontrada en = ", x_bar)
else
	println("Orbita de periodo 3 no encontrada")
end
println((isguaranteed(Y), isguaranteed(Z₁), isguaranteed(Z₂)))

println("")
println("***************************ESTABILIDAD DE LA ORBITA***************************")

# Definimos G = f o f o f, y DG = D(f o f o f) = Df((fof)(x0))Df(f(x0))Df(x0) = Df(x2)Df(x1)Df(x0)
G(x) = f(f(f(x,μ),μ),μ)
DG(x) = Df(f(f(x,μ),μ),μ) * Df(f(x,μ),μ) * Df(x,μ)
println(estabilidad(DG(ix_bar[1])))
