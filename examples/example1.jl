using IntervalArithmetic

a = 0.5
b = 0.1

@show(a)
@show(b)
@show(interval(a))
@show(interval(b))
@show(interval(a) + interval(b))
@show(interval(a) - interval(b))
@show(interval(a) * interval(b))
@show(interval(a) / interval(b))