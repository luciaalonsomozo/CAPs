using RadiiPolynomial
using CSV
using DataFrames

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
            initial_guess = rand(n)
            println("New initial guess = ", initial_guess)
        end
    end

    return 0, false
end

function agregar_fila_csv(archivo::String, mu::String, n::Int, texto::String)
    df = DataFrame(nu = [mu], periodo = [n], estabilidad = [texto])
    CSV.write(archivo, df; append=true, writeheader=!isfile(archivo))
end

function run_main(iμ_str::String, n::Int)
    imu = eval(Meta.parse("@interval($iμ_str)"))
    mu = mid(imu)

    F_DF(x) = (F(x, mu, n), DF(x, mu, n))
    x_bar, success = newton_orbit_algorithm(F_DF, n; error_threshold=ERROR, max_iter=20)

    if !success
        println("No se encontró órbita de periodo $n para μ = $mu")
        return
    end

    println("x_bar = ", x_bar)

    A = inv(DF(x_bar, mu, n))
    ix_bar = interval(x_bar)
    iA = interval(A)

    Y = norm(iA * F(ix_bar, imu, n), 1)
    Z1 = opnorm(interval(I) - iA * DF(ix_bar, imu, n), 1)
    R = Inf
    Z2 = interval(2) * imu * opnorm(iA, 1)
    ie = interval_of_existence(Y, Z1, Z2, R)
    r = inf(ie)

    println("r = ", r)
    println((isguaranteed(Y), isguaranteed(Z1), isguaranteed(Z2)))

    enclosures = [interval(ix_bar[i], r; format = :midpoint) for i in eachindex(ix_bar)]
    if isdisjoint_interval(enclosures...)
        println("Órbita de periodo $n encontrada en = ", ix_bar)
    else
        println("Órbita de periodo $n no encontrada porque los intervalos no son disjuntos")
        return
    end

    println("***************************ESTABILIDAD DE LA ÓRBITA***************************")
    agregar_fila_csv("datos_logistica_4.csv", iμ_str, n, estabilidad(DG(ix_bar, imu)))
end

function main()
    for mu_str in string.(3.4:0.0005:4.0)
        for n in 2:11
            for _ in 1:5  # 5 intentos por combinación
                println("Running for μ = $mu_str, n = $n")
                run_main(mu_str, n)
            end
        end
    end

end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end