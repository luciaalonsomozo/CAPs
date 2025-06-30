using Plots
using RadiiPolynomial
using CSV
using DataFrames

include("estabilidad.jl")
using .Estabilidad

const ERROR = 1.0e-3

f(x, b, k) = b + x - k / (ExactReal(1) + exp(x))
Df(x, k) = ExactReal(1) + k * exp(x) / (exp(x) + ExactReal(1))^2

function F(x, b, k, n)
    v = zeros(eltype(x), n)
    for i in 1:n
        v[i] = f(x[i], b, k) - x[mod1(i+1, n)]
    end
    return v
end

function DF(x, k, n)
    M = zeros(eltype(x), n, n)
    for i in 1:n
        M[i,i] = Df(x[i], k)
        M[i,mod1(i+1, n)] = ExactReal(-1)
    end
    return M
end

function DG(x, k)
    result = ExactReal(1)
    for i in eachindex(x)
        result *= Df(x[i], k)
    end
    return result
end

function newton_orbit_algorithm(F_DF, n; error_threshold=1e-6, max_iter=100)
    iter = 0
    initial_guess = 10 * rand(n)
    println("initial_guess = ", initial_guess)
    x_bar = initial_guess
    success = false

    while !success && iter < max_iter
        iter += 1
        try
            x_bar, success = newton(F_DF, initial_guess, maxiter=20)
        catch
        end

        max_diff = maximum(x_bar) - minimum(x_bar)
        if max_diff > error_threshold && success
            println("Not a fixed point.")
            return x_bar, true
        else
            println(success ? "Fixed point found: $x_bar" : "Failed iteration.")
            initial_guess = 10 * rand(n)
            println("initial_guess = ", initial_guess)
        end
    end

    return 0, false
end

function agregar_fila_csv(archivo::String, a::String, b::String, n::Int, texto::String)
    df = DataFrame(b = [a], k = [b], periodo = [n], estabilidad = [texto])
    # Usamos `append=true` para no sobrescribir el archivo
    CSV.write(archivo, df; append=true, writeheader=!isfile(archivo))
end

function run_main(b_str::String, k_str::String, n::Int)
    ib = eval(Meta.parse("@interval($b_str)"))
    ik = eval(Meta.parse("@interval($k_str)"))

    b = mid(ib)
    k = mid(ik)
    F_DF(x) = (F(x, b, k, n), DF(x, k, n))
    x_bar, success = newton_orbit_algorithm(F_DF, n; error_threshold=ERROR, max_iter=10)

    if !success
        println("No se encontró una órbita de periodo $n para b = $b y k = $k.")
        return
    end

    println("x_bar = ", x_bar)

    A = inv(DF(x_bar, k, n))
    ik = interval(k)
    ib = interval(b)
    ix_bar = interval(x_bar)
    iA = interval(A)

    Y = norm(iA * F(ix_bar, ib, ik, n), 1)
    Z₁ = opnorm(interval(I) - iA * DF(ix_bar, ik, n), 1)
    R = Inf
    ik_modulo = -interval(ik)
    L = sqrt(interval(3)) / interval(18) * ik_modulo
    Z₂ = L * opnorm(iA, 1)

    ie = interval_of_existence(Y, Z₁, Z₂, R)
    r = inf(ie)

    println("r = ", r)
    println((isguaranteed(Y), isguaranteed(Z₁), isguaranteed(Z₂)))

    enclosures = [interval(ix_bar[i], r; format=:midpoint) for i in eachindex(ix_bar)]

    if isdisjoint_interval(enclosures...)
        println("Órbita de periodo $n encontrada en = ", ix_bar)
    else
        println("Órbita de periodo $n no encontrada porque los intervalos no son disjuntos")
        return
    end

    println("***************************ESTABILIDAD DE LA ORBITA***************************")
    agregar_fila_csv("datos_dp_11_06_aumentando_periodo.csv", b_str, k_str, n, estabilidad(DG(ix_bar, ik)))
end

function main()
    for b_val in -20.0:0.25:-2.0
        for k_val in -5.0:-0.25:-45.0
            for n in 2:10
                println("Running for b = $b_val, k = $k_val, n = $n")
                run_main(string(b_val), string(k_val), n)
            end
        end
    end

end