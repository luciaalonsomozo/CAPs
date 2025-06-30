# orbita_general.jl

using Plots
using RadiiPolynomial

include("estabilidad.jl")  # Incluir el archivo
using .Estabilidad  # Importar el módulo

include("graficar_orbita_logistica.jl")
using .Graficar

ERROR = 1.0e-3

f(x, μ) = μ * x * (ExactReal(1) - x) # logistic map
Df(x, μ) = μ * (ExactReal(1) - ExactReal(2)*x) # derivative of the logistic map w.r.t. x

function F(x, μ)
    v = zeros(eltype(x), n)
    for i in 1:n
        v[i] = f(x[i], μ) - x[mod1(i+1, n)]
    end
	return v
end

function DF(x, μ) # derivative of the zero-finding problem w.r.t. x
	M = zeros(eltype(x), n, n)
    for i in 1:n
        M[i,i] = Df(x[i], μ)
        M[i,mod1(i+1, n)] = ExactReal(-1)
    end
	return M
end

# Función para verificar si todos los intervalos son disjuntos (pares)
function are_disjoint(enclosures)
    for i in eachindex(enclosures)
        for j in (i+1):length(enclosures)
            if !isdisjoint_interval(enclosures[i], enclosures[j])
                return false
            end
        end
    end
    return true
end

function DG(x, μ)
    result = ExactReal(1)
    for i in eachindex(x)
        println("Para x = ", x[i], " ,tenemos = ", Df(x[i], μ))
        result = result * Df(x[i], μ)
    end
    return result
end

function newton_orbit_algorithm(F_DF, error_threshold=1e-6, max_iter=100)
    # Start the Newton method with an initial guess
    newton_orbit = false
    iter = 0

    initial_guess = rand(n)
    println("initial_guess = ", initial_guess)

    while !newton_orbit && iter < max_iter
        
        iter += 1

        # Apply Newton's method
        x_bar, success = newton(F_DF, initial_guess, maxiter = 50)

        # Check if the components of x_bar are within the error threshold
        # Assuming x_bar is a vector and we want the components to be small enough
        max_diff = maximum(x_bar) - minimum(x_bar)  # Find the difference between max and min component
        # isdisjoint_interval(interval(x_bar)...)
        if max_diff > error_threshold && success
            newton_orbit = true  # Converged if max_error is below threshold
            println("Not a fixed point.")
            return x_bar, true
        elseif !success
            # Update initial_guess for the next iteration
            initial_guess = rand(n)
            println("Values got: ", x_bar)
            println("Continuing to next iteration with new initial_guess = ", initial_guess)
        else
            println("Fixed point found: ", x_bar)
            # Update initial_guess for the next iteration
            initial_guess = rand(n)
            println("Continuing to next iteration with new initial_guess = ", initial_guess)
        end
    end

    if !newton_orbit
        return 0, false
    end
end


if length(ARGS) < 2
    println("Error: debes pasar al menos dos argumentos")
    exit(1)  # Termina la ejecución con error
end

n = parse(Int, ARGS[1])
iμ = parse(Interval{Float64}, ARGS[2])
#@show in_interval(39//10, iμ)
#= This returns false.
μ = parse(Float64, ARGS[2])
a = interval(μ)
@show in_interval(39//10, a)
=#
μ = mid(iμ)

F_DF(x) = (F(x, μ), DF(x, μ)) # returns `F` and `DF' at `μ`
x_bar, success = newton_orbit_algorithm(F_DF, ERROR, 5)

if success == false
    println("No se encontró una órbita de periodo $n en el máximo número de iteraciones.")
    exit(1)
end

println("x_bar = ", x_bar)

A = inv(DF(x_bar, μ))
ix_bar = interval(x_bar)
iA = interval(A)
Y = norm(iA * F(ix_bar, iμ), 1)
Z₁ = opnorm(interval(I) - iA * DF(ix_bar, iμ), 1)
R = Inf
Z₂ = interval(2) * iμ * opnorm(iA, 1)
ie = interval_of_existence(Y, Z₁, Z₂, R)
r = inf(ie)

println("r = ", r)

# Crear automáticamente los intervalos:
enclosures = [interval(ix_bar[i], r; format = :midpoint) for i in eachindex(ix_bar)]

# Comprobación de disjunción:
if isdisjoint_interval(enclosures...)#are_disjoint(enclosures)
    println("Órbita de periodo $n encontrada en = ", ix_bar)
else
    println("Órbita de periodo $n no encontrada")
end

println((isguaranteed(Y), isguaranteed(Z₁), isguaranteed(Z₂)))

println("")
println("***************************ESTABILIDAD DE LA ORBITA***************************")

println(estabilidad(DG(enclosures, iμ)))

# graficar(f, x_bar[1], μ)