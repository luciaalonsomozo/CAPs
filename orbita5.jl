# orbita5.jl

using Plots
using RadiiPolynomial

include("estabilidad.jl")  # Incluir el archivo
using .Estabilidad  # Importar el módulo

f(x, μ) = μ * x * (1 - x) # logistic map

function F(x, μ)
    v = zeros(eltype(x), 5)
    v[1] = f(x[1], μ) - x[2]
	v[2] = f(x[2], μ) - x[3]
	v[3] = f(x[3], μ) - x[4]
    v[4] = f(x[4], μ) - x[5]
    v[5] = f(x[5], μ) - x[1]
	return v
end

Df(x, μ) = μ * (1 - 2x) # derivative of the logistic map w.r.t. x

function DF(x, μ) # derivative of the zero-finding problem w.r.t. x
	M = zeros(eltype(x), 5, 5)
	M[1,1] = Df(x[1], μ)
	M[1,2] = -1
	M[2,2] = Df(x[2], μ)
	M[2,3] = -1
	M[3,3] = Df(x[3], μ)
    M[3,4] = -1
    M[4,4] = Df(x[4], μ)
    M[4,5] = -1
    M[5,5] = Df(x[5], μ)
    M[5,1] = -1
	return M
end

initial_guess = [0.24, 0.682, 0.81, 0.57, 0.916]  # Estimación numérica
μ  = 3.74 # ponemos este nu para buscar órbita de periodo 5
F_DF(x) = (F(x, μ), DF(x, μ)) # returns `F` and `DF` at `μ = 3.74`
x_bar, success = newton(F_DF, initial_guess; tol=1.0e-12, maxiter = 50, verbose=true)
println(x_bar)


A = inv(DF(x_bar, μ))
# since 374/100 is not representable as a floating-point number
iμ = interval(374)/interval(100)
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
enclosure_x_tilde_4 = interval(ix_bar[4], r; format = :midpoint)
enclosure_x_tilde_5 = interval(ix_bar[5], r; format = :midpoint)

if isdisjoint_interval(enclosure_x_tilde_1, enclosure_x_tilde_2, enclosure_x_tilde_3, enclosure_x_tilde_4, enclosure_x_tilde_5)
	println("Orbita de periodo 5 encontrada en = ", x_bar)
else
	println("Orbita de periodo 5 no encontrada")
end

println((isguaranteed(Y), isguaranteed(Z₁), isguaranteed(Z₂)))

println("")
println("***************************ESTABILIDAD DE LA ORBITA***************************")

# Definimos G = f o f o f o f o f
function DG(x)
	return Df(x[1],μ)*Df(x[2],μ)*Df(x[3],μ)*Df(x[4],μ)*Df(x[5],μ)
end

println(estabilidad(DG(ix_bar)))



