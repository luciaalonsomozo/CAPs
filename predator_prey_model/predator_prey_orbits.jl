using Plots
using RadiiPolynomial, Random

include("../stability/stability.jl")
using .Stability

const ERROR = 1.0e-3

# Define the function and its derivative
f(x, b, k) = b + x - k / (ExactReal(1) + exp(x))
Df(x, k) = ExactReal(1) + k * exp(x) / (exp(x) + ExactReal(1))^2

# Defines the vector field F corresponding to the orbit
function F(x, b, k, n)
    v = zeros(eltype(x), n)
    for i in 1:n
        v[i] = f(x[i], b, k) - x[mod1(i+1, n)]
    end
    return v
end

# Computes the Jacobian matrix DF of F
function DF(x, k, n)
    M = zeros(eltype(x), n, n)
    for i in 1:n
        M[i,i] = Df(x[i], k)
        M[i,mod1(i+1, n)] = ExactReal(-1)
    end
    return M
end

# Computes the product of derivatives along the orbit
function DG(x, k)
    result = ExactReal(1)
    for i in eachindex(x)
        result *= Df(x[i], k)
    end
    return result
end

# Algorithm to find a periodic orbit using Newton’s method
function newton_orbit_algorithm(F_DF, n; error_threshold=1e-6, max_iter=100)
    iter = 0
    initial_guess = 10 * rand(n)
    x_bar = initial_guess
    success = false

    while !success && iter < max_iter
        iter += 1
        try
            x_bar, success = newton(F_DF, initial_guess, maxiter=50)
        catch
        end

        max_diff = maximum(x_bar) - minimum(x_bar)
        if max_diff > error_threshold && success
            println("Not a fixed point.")
            return x_bar, true
        else
            println(success ? "Fixed point found: $x_bar" : "Failed iteration.")
            initial_guess = 10 * rand(n)
        end
    end

    return 0, false
end

# Main function that validates the orbit and checks its stability
function main(b_str::String, k_str::String, n::Int)
    ib = eval(Meta.parse("@interval($b_str)"))
    ik = eval(Meta.parse("@interval($k_str)"))

    b = mid(ib)
    k = mid(ik)

    F_DF(x) = (F(x, b, k, n), DF(x, k, n))
    x_bar, success = newton_orbit_algorithm(F_DF, n; error_threshold=ERROR, max_iter=1000)

    if !success
        println("No periodic orbit of period $n was found.")
        return
    end

    println("x_bar = ", x_bar)

    A = inv(DF(x_bar, k, n))
    ix_bar = interval(x_bar)
    iA = interval(A)

    Y = norm(iA * F(ix_bar, ib, ik, n), 1)
    Z₁ = opnorm(interval(I) - iA * DF(ix_bar, ik, n), 1)
    R = Inf
    ik_modulo = -interval(ik)
    L = interval(2) * sqrt(interval(3)) / interval(18) * ik_modulo
    Z₂ = L * opnorm(iA, 1)

    println("Y = ", Y)
    println("Z1 = ", Z₁)
    println("Z2 = ", Z₂)

    ie = interval_of_existence(Y, Z₁, Z₂, R)
    r = inf(ie)

    println("r = ", r)
    println((isguaranteed(Y), isguaranteed(Z₁), isguaranteed(Z₂)))

    enclosures = [interval(ix_bar[i], r; format=:midpoint) for i in eachindex(ix_bar)]

    if isdisjoint_interval(enclosures...)
        println("Periodic orbit of period $n found in = ", ix_bar)
    else
        println("Periodic orbit of period $n not found because the intervals are not disjoint.")
        return
    end

    println("***************************STABILITY OF THE ORBIT***************************")
    println(stability(DG(enclosures, ik)))
    println(DG(enclosures, ik))
end
