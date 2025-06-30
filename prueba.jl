using IntervalArithmetic

iμ = parse(Interval{Float64}, ARGS[1])
@show in_interval(39//10, iμ)

μ = parse(Float64, ARGS[1])
a = interval(μ)
@show in_interval(39//10, a)

mu_val = 3.9
mu = ExactReal(3.9)
imu = interval(mu)
@show in_interval(39//10, imu)