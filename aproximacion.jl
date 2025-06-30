# aproximacion.jl

using Plots
using RadiiPolynomial
using LinearAlgebra


f1(x) = exp(x) # ecuación diferencial 1
f2(x) = exp(x) - exp(-10)# ecuación diferencial 2

Df1(x) = exp(x)
Df2(x) = exp(x)

# ecuación 2

println("Ecuación f2(x) = exp(x) - exp(-10) =====> debería dar cierto")

initial_guess = 0

f_df_2(x) = (f2(x), Df2(x))
x_bar2, success2 = newton(f_df_2, initial_guess, maxiter = 100)
println(x_bar2)

A = inv(Df2(x_bar2))
ix_bar2 = interval(x_bar2)
iA = interval(A)
Y = norm(iA * f2(ix_bar2), 1)
Z₁ = opnorm(I - iA * Df2(ix_bar2), 1)
R = ExactReal(5)
Z₂ = interval(exp(R))
ie = interval_of_existence(Y, Z₁, Z₂, R)
r = inf(ie)

println("r = ", r)

println((isguaranteed(Y), isguaranteed(Z₁), isguaranteed(Z₂)))

# ecuación 1

println("Ecuación f1(x) = exp(x) =====> debería dar false")


initial_guess = 0.1
f_df_1(x) = (f1(x), Df1(x))
x_bar1, success1 = newton(f_df_1, initial_guess, maxiter = 100)
println(x_bar1)

A = inv(Df1(x_bar1))
ix_bar1 = interval(x_bar1)
iA = interval(A)
Y = norm(iA * f1(ix_bar1), 1)
Z₁ = opnorm(I - iA * Df1(ix_bar1), 1)
R = ExactReal(5)
Z₂ = interval(exp(R))
ie = interval_of_existence(Y, Z₁, Z₂, R)
r = inf(ie)

println("r = ", r)

println((isguaranteed(Y), isguaranteed(Z₁), isguaranteed(Z₂)))