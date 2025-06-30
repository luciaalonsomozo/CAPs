using RadiiPolynomial

f(x, μ) = μ * x * (ExactReal(1) - x)
Df(x, μ) = μ * (ExactReal(1) - ExactReal(2)*x)

function F(x, μ, space)
   return project(f(x, μ) - x, space)
end

function DF(x, μ, domain, codomain)
    return project(Multiplication(Df(x, μ) - ExactReal(1)), domain, codomain)
end

function main(N = 5)
    μ = Sequence(Taylor(1), [3.9, 0.1])
    iμ = Sequence(Taylor(1), [I"3.9", I"0.1"])

    F_DF(x) = (F(x, μ, space(x)), DF(x, μ, space(x), space(x)))

    initial_guess = ones(ComplexF64, Taylor(N)) # do not start at zero, it's a trivial fixed-point
    x_bar, success = newton(F_DF, initial_guess)

    println("x_bar = ", x_bar)

    # A_k = inv(DF(x_bar, μ, space(x_bar), space(x_bar)))
    A_k = inv(- μ * x_bar)
    ix_bar = interval(x_bar)
    iA_k = interval(A_k)

    Y = norm(iA_k * F(ix_bar, iμ, Taylor(2order(x_bar)+order(μ))), 1)
    # IMPORTANT to get the proper range of F : 2order(x_bar)+order(μ), the 2 comes from the quadratic nonlinearity in x
    Z₁ = opnorm(interval(I) - iA_k * DF(ix_bar, iμ, space(ix_bar), space(ix_bar)), 1)
    R = Inf
    Z₂ = interval(2) * norm(iμ, 1) * opnorm(iA_k, 1)
    ie = interval_of_existence(Y, Z₁, Z₂, R)
    r = inf(ie)

    println("r = ", r)

    println((isguaranteed(Y), isguaranteed(Z₁), isguaranteed(Z₂)))
end

main()
