using RadiiPolynomial

include("estabilidad.jl")
using .Estabilidad

const ERROR = 1.0e-3

f(x, μ) = μ * x * (ExactReal(1) - x)
Df(x, μ) = μ * (ExactReal(1) - ExactReal(2)*x)

function F(x, μ, n)
    v = zeros(eltype(x), n)
    for i in 1:n
        v[i] = f(x[i], μ) - x[mod1(i+1, n)]
    end
    return v
end

function DF(x, μ, n)
    M = zeros(eltype(x), n, n)
    for i in 1:n
        M[i,i] = Df(x[i], μ)
        M[i,mod1(i+1, n)] = ExactReal(-1)
    end
    return M
end

function DG(x, μ)
    result = ExactReal(1)
    for i in eachindex(x)
        result *= Df(x[i], μ)
    end
    return result
end

function newton_orbit_algorithm(F_DF, n; error_threshold=1e-6, max_iter=100)
    iter = 0
    initial_guess = rand(n)
    println("initial_guess = ", initial_guess)

    while iter < max_iter
        iter += 1
        x_bar, success = newton(F_DF, initial_guess, maxiter = 50)
        max_diff = maximum(x_bar) - minimum(x_bar)
        if max_diff > error_threshold && success
            println("Not a fixed point.")
            return x_bar, true
        else
            println(success ? "Fixed point found: $x_bar" : "Failed iteration.")
            initial_guess = rand(n)
            println("New initial guess = ", initial_guess)
        end
    end
    return 0, false
end

function main(iμ_str::String, n::Int)
    iμ = eval(Meta.parse("@interval($iμ_str)"))
    μ = mid(iμ)

    F_DF(x) = (F(x, μ, n), DF(x, μ, n))
    x_bar, success = newton_orbit_algorithm(F_DF, n; error_threshold=ERROR, max_iter=1000)

    if !success
        println("No se encontró una órbita de periodo $n.")
        return
    end

    println("x_bar = ", x_bar)

    A = inv(DF(x_bar, μ, n))
    ix_bar = interval(x_bar)
    iA = interval(A)
    Y = norm(iA * F(ix_bar, iμ, n), 1)
    Z₁ = opnorm(interval(I) - iA * DF(ix_bar, iμ, n), 1)
    R = Inf
    Z₂ = interval(2) * iμ * opnorm(iA, 1)
    ie = interval_of_existence(Y, Z₁, Z₂, R)
    r = inf(ie)

    println("r = ", r)

    enclosures = [interval(ix_bar[i], r; format = :midpoint) for i in eachindex(ix_bar)]
    if isdisjoint_interval(enclosures...)
        println("Órbita de periodo $n encontrada en = ", ix_bar)
    else
        println("Órbita de periodo $n no encontrada")
        return
    end

    println((isguaranteed(Y), isguaranteed(Z₁), isguaranteed(Z₂)))
    println("***************************ESTABILIDAD DE LA ORBITA***************************")
    println(estabilidad(DG(enclosures, iμ)))
end
