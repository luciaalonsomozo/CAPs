# orbita2.jl

using Plots
using RadiiPolynomial
using LinearAlgebra

include("estabilidad.jl")  # Incluir el archivo
using .Estabilidad  # Importar el módulo

f(x, b, k) = b + x - k/(1 + exp(x)) # ecuación diferencial

function F(x, b, k)
    v = zeros(eltype(x), 2)
    v[1] = f(x[1], b, k) - x[2]
	v[2] = f(x[2], b, k) - x[1]
	return v
end

# Podemos desprendernos de b para la derivada
Df(x,k) = 1 + k*exp(x)/(exp(x)+1)^2 # derivative of the map 

function DF(x, k) # derivative of the zero-finding problem w.r.t. x
	M = zeros(eltype(x), 2, 2)
	M[1,1] = Df(x[1], k)
	M[1,2] = -1
    M[2,1] = -1
	M[2,2] = Df(x[2], k)
	return M
end

println("")
println("***************************ORBITA***************************")
# con este encuentra la orbita de periodo 2
initial_guess = [0.1, 5.0]
# initial_guess = [5.0, 6.0]
b = -3
k = -10
F_DF(x) = (F(x, b, k), DF(x, k)) # returns `F` and `DF` at `b = -10` and `k = -3`
x_bar, success = newton(F_DF, initial_guess, maxiter = 100)

A = inv(DF(x_bar, k))
ik = interval(-10)
ib = interval(-3)
ix_bar = interval(x_bar)
iA = interval(A)
println("ik = ", ik)
println("ib = ", ib)
println("ix_bar = ", ix_bar)
println("x_bar = ", x_bar)
println("A = ", A)
println("iA = ", iA)
Y = norm(iA * F(ix_bar, ib, ik), 1)
Z₁ = opnorm(I - iA * DF(ix_bar, ik), 1)
R = Inf
L = interval(sqrt(3))/interval(18)
Z₂ = interval(2) * L * interval(10) * opnorm(iA, 1)
println("Y = ", Y)
println("Z₁ = ", Z₁)
println("Z₂ = ", Z₂)
ie = interval_of_existence(Y, Z₁, Z₂, R)
r = inf(ie)

println("r = ", r)

enclosure_x_tilde_1 = interval(ix_bar[1], r; format = :midpoint)
enclosure_x_tilde_2 = interval(ix_bar[2], r; format = :midpoint)

if isdisjoint_interval(enclosure_x_tilde_1, enclosure_x_tilde_2)
	println("Orbita de periodo 2 encontrada en = ", x_bar)
else
	println("Orbita de periodo 2 no encontrada")
end
println((isguaranteed(Y), isguaranteed(Z₁), isguaranteed(Z₂)))

println("")
println("***************************ESTABILIDAD DE LA ORBITA***************************")

# Definimos G = f o f, y DG = D(f o f) = Df(f(x0))Df(x0) = Df(x1)Df(x0)
G(x) = f(f(x, b, k), b, k)
DG(x) = Df(f(x, b, k), k) * Df(x, k)
println(estabilidad(DG(ix_bar[1])))

println("")
println("***************************PUNTO FIJO***************************")
# con este encuentra el punto fijo
initial_guess = [0.5, 1.0]
x_puntofijo, success = newton(F_DF, initial_guess, maxiter = 100)
ix_puntofijo = interval(x_puntofijo[1])
println("Punto fijo encontrado en = ", ix_puntofijo)

println("")
println("***************************ESTABILIDAD DEL PUNTO FIJO***************************")
println(estabilidad(Df(ix_puntofijo, k)))