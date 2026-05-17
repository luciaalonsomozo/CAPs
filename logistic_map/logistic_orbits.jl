using RadiiPolynomial, Random
using CSV
using DataFrames
using LinearAlgebra

include("../stability/stability.jl")
using .Stability

const ERROR = 1.0e-6

# Plain Float64 versions for Newton iteration
f(x::Float64, μ::Float64)  = μ * x * (1.0 - x)
Df(x::Float64, μ::Float64) = μ * (1.0 - 2.0*x)

# Generic versions that work with intervals AND Float64
f_generic(x, μ)  = μ * x * (1 - x)
Df_generic(x, μ) = μ * (1 - 2*x)

function F(x, μ, n)
    v = zeros(eltype(x), n)
    for i in 1:n
        v[i] = f_generic(x[i], μ) - x[mod1(i+1, n)]
    end
    return v
end

function DF(x, μ, n)
    M = zeros(eltype(x), n, n)
    for i in 1:n
        M[i, i]            = Df_generic(x[i], μ)
        M[i, mod1(i+1, n)] = -one(eltype(x))
    end
    return M
end

function DG(x, μ)
    result = one(eltype(x))
    for i in eachindex(x)
        result *= Df_generic(x[i], μ)
    end
    return result
end

function estabilidad(landa)
    if isstrictless(interval(1), abs(landa))
        return "U"
    elseif in_interval(1, abs(landa))
        return "N"
    else
        return "S"
    end
end

function newton_orbit_algorithm(mu::Float64, n::Int; error_threshold=1e-6, max_iter=100)
    iter    = 0
    success = false
    x_bar   = rand(n)

    while !success && iter < max_iter
        iter += 1
        println("initial_guess = ", x_bar)

        Fval = [f(x_bar[i], mu) - x_bar[mod1(i+1, n)] for i in 1:n]
        J    = zeros(Float64, n, n)
        for i in 1:n
            J[i, i]            = Df(x_bar[i], mu)
            J[i, mod1(i+1, n)] = -1.0
        end

        local x_new
        try
            x_new = x_bar - J \ Fval
        catch
            println("Linear solve failed, retrying.")
            x_bar = rand(n)
            continue
        end

        for _ in 1:200
            Fval = [f(x_new[i], mu) - x_new[mod1(i+1, n)] for i in 1:n]
            res  = maximum(abs, Fval)
            if res < 1e-12
                success = true
                break
            end
            J = zeros(Float64, n, n)
            for i in 1:n
                J[i, i]            = Df(x_new[i], mu)
                J[i, mod1(i+1, n)] = -1.0
            end
            try
                x_new = x_new - J \ Fval
            catch
                break
            end
        end

        if success
            max_diff = maximum(x_new) - minimum(x_new)
            if max_diff > error_threshold
                println("Orbit candidate found: ", x_new)
                return x_new, true
            else
                println("Collapsed to fixed point, retrying.")
                success = false
            end
        else
            println("Failed iteration.")
        end

        x_bar = rand(n)
    end

    return zeros(n), false
end

function agregar_fila_csv(archivo::String, mu::String, n::Int, texto::String)
    df = DataFrame(nu=[mu], periodo=[n], estabilidad=[texto])
    CSV.write(archivo, df; append=true, writeheader=!isfile(archivo))
end

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

# Lift a Float64 to a thin guaranteed interval using prevfloat/nextfloat for rigorous bounds
float_to_interval(x::Float64) = interval(prevfloat(x), nextfloat(x))

function run_main(iμ_str::String, n::Int)
    imu = eval(Meta.parse("@interval($iμ_str)"))
    mu  = mid(imu)

    x_bar, success = newton_orbit_algorithm(mu, n; error_threshold=ERROR, max_iter=100)

    if !success
        println("No periodic orbit of period $n found for μ = $mu")
        return false
    end

    println("x_bar = ", x_bar)
    println("***************************STABILITY OF THE ORBIT***************************")

    # Lift Float64 orbit points to guaranteed intervals using directed rounding
    ix_bar = [float_to_interval(x_bar[i]) for i in eachindex(x_bar)]

    # Compute A = inv(DF) in Float64, then lift to guaranteed intervals
    A  = inv(DF(x_bar, mu, n))
    iA = [float_to_interval(A[i,j]) for i in 1:n, j in 1:n]

    Y  = norm(iA * F(ix_bar, imu, n), 1)
    Z1 = opnorm(interval.(I(n)) - iA * DF(ix_bar, imu, n), 1)
    R  = Inf
    Z2 = interval(2) * imu * opnorm(iA, 1)

    println("Y  = ", Y)
    println("Z1 = ", Z1)
    println("Z2 = ", Z2)

    ie, success_ie = interval_of_existence(Y, Z1, Z2, R)
    if !success_ie
        println("interval_of_existence failed.")
        return false
    end

    r = inf(ie)
    println("r  = ", r)

    enclosures = [interval(ix_bar[i], r; format=:midpoint) for i in eachindex(ix_bar)]

    if are_disjoint(enclosures)
        λ    = DG(ix_bar, imu)
        stab = estabilidad(λ)
        println("Period-$n orbit found at = ", ix_bar, "  multiplier = ", λ)
        agregar_fila_csv("80_logistic_map.csv", iμ_str, n, stab)
        return true
    else
        println("Period-$n orbit not verified — enclosures not disjoint.")
        return false
    end
end

# Main function reads mu_str and n parameters directly from terminal inputs
function main(mu_str::String, n::Int)
    for attempt in 1:50
        if run_main(mu_str, n) == true
            println("Successfully verified orbit on attempt $attempt.")
            break
        end
    end
end

# ==============================================================================
# SCAN LOGIC FOR GENERATING THOUSANDS OF ORBITS (COMMENTED OUT)
# ==============================================================================
# function main()
#     # Comprehensive sequence checking period-doubling cascade spikes first
#     periods = [2, 4, 8, 16, 3, 6, 5, 7, 9, 10, 11, 12, 13, 14, 15,
#                17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
#                31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44,
#                45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58,
#                59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72,
#                73, 74, 75, 76, 77, 78, 79, 80, 81]
#
#     # Parametric sweep over the interesting chaotic window of the Logistic Map
#     for mu_val in 3.4:0.0005:4.0
#         mu_str = string(round(mu_val, digits=4))
#         
#         for n in periods
#             for attempt in 1:5   # Up to 5 verification attempts per (μ, n) configuration
#                 if run_main(mu_str, n) == true
#                     break        # Orbit verified! Move forward to next period size
#                 end
#             end
#         end
#     end
# end
# ==============================================================================