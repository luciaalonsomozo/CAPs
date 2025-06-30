# orbita_general.jl

using Plots
using RadiiPolynomial, Random

include("estabilidad.jl")  # Incluir el archivo
using .Estabilidad  # Importar el módulo

include("graficar_orbita_depr_presa.jl")
using .Graficar

rng = MersenneTwister()


if length(ARGS) < 3
    println("Error: debes pasar al menos tres argumentos")
    exit(1)  # Termina la ejecución con error
end

n = parse(Int, ARGS[1]) # Tamaño de la órbita (2^n)
n = 2^n
b = parse(Float64, ARGS[2]) # Parametro b
k = parse(Float64, ARGS[3]) # Parametro k

ERROR = 1.0e-3

f(x, b, k) = b + x - k/(ExactReal(1) + exp(x)) # ecuación diferencial
Df(x,k) = ExactReal(1) + k*exp(x)/(exp(x)+ExactReal(1))^2 # derivative of the map 

function F(x, b, k)
    v = zeros(eltype(x), n)
    for i in 1:n
        v[i] = f(x[i], b, k) - x[mod1(i+1, n)]
    end
	return v
end

function DF(x, k) # derivative of the zero-finding problem w.r.t. x
	M = zeros(eltype(x), n, n)
    for i in 1:n
        M[i,i] = Df(x[i], k)
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

function DG(x, k)
    result = ExactReal(1)
    for i in eachindex(x)
        println("Para x = ", x[i], " ,tenemos = ", Df(x[i], k))
        result = result * Df(x[i], k)
    end
    return result
end

function newton_orbit_algorithm(F_DF, error_threshold=1e-6, max_iter=100)
    # Start the Newton method with an initial guess
    newton_orbit = false
    iter = 0

    initial_guess = 10* rand(n)
    # initial_guess = [0.1, 5.0]
    println("initial_guess = ", initial_guess)

    x_bar = initial_guess
    success = false

    while !newton_orbit && iter < max_iter
        
        iter += 1

        # Apply Newton's method
        try
            x_bar, success = newton(F_DF, initial_guess, maxiter = 50)
        catch
        end


        # Check if the components of x_bar are within the error threshold
        # Assuming x_bar is a vector and we want the components to be small enough
        max_diff = maximum(x_bar) - minimum(x_bar)  # Find the difference between max and min component
        if max_diff > error_threshold && success
            newton_orbit = true  # Converged if max_error is below threshold
            println("Not a fixed point.")
            return x_bar, true
        elseif !success
            # Update initial_guess for the next iteration
            @show initial_guess = 10*rand(n)
            println("Values got: ", x_bar)
            println("Continuing to next iteration with new initial_guess = ", initial_guess)
        else
            println("Fixed point found: ", x_bar)
            # Update initial_guess for the next iteration
            @show initial_guess = 10*rand(n)
            println("Continuing to next iteration with new initial_guess = ", initial_guess)
        end
    end

    if !newton_orbit
        return 0, false
    end
end

F_DF(x) = (F(x, b, k), DF(x, k)) # returns `F` and `DF' at `μ`
x_bar, success = newton_orbit_algorithm(F_DF, ERROR, 50)

if success == false
    println("No se encontró una órbita de periodo $n en el máximo número de iteraciones.")
    exit(1)
end

println("x_bar = ", x_bar)

A = inv(DF(x_bar, k))
ik = interval(k)
ib = interval(b)
ix_bar = interval(x_bar)
iA = interval(A)
Y = norm(iA * F(ix_bar, ib, ik), 1)
Z₁ = opnorm(interval(I)- iA * DF(ix_bar, ik), 1)
R = Inf
ik_modulo = -interval(ik)
L = interval(sqrt(3))/interval(18) * ik_modulo
# Z₂ = interval(2) * L * interval(10) * opnorm(iA, 1)
Z₂ = L * opnorm(iA, 1)

ie = interval_of_existence(Y, Z₁, Z₂, R)
r = inf(ie)

println("r = ", r)

# Crear automáticamente los intervalos:
enclosures = [interval(ix_bar[i], r; format = :midpoint) for i in eachindex(ix_bar)]

# Comprobación de disjunción:
if are_disjoint(enclosures)
    println("Órbita de periodo $n encontrada en = ", ix_bar)
else
    println("Órbita de periodo $n no encontrada porque los intervalos no son disjuntos")
end

println((isguaranteed(Y), isguaranteed(Z₁), isguaranteed(Z₂)))

println("")
println("***************************ESTABILIDAD DE LA ORBITA***************************")

println(estabilidad(DG(ix_bar, ik)))

# graficar(f, x_bar[1], b, k)

