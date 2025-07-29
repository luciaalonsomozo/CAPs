using IntervalArithmetic

a = 0.5
b = 0.1
c = interval(-0.5, 1.5)

@show(interval(a))
@show(interval(b))
@show(interval(c))
@show(isequal_interval(interval(c)*(interval(a)+interval(b)), interval(c)*interval(a) + interval(c)*interval(b)))
