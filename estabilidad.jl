module Estabilidad

export estabilidad

using LinearAlgebra
using IntervalArithmetic
using RadiiPolynomial

function estabilidad(landa)

    if isstrictless(interval(1), abs(landa))
        return "I"
    elseif in_interval(1, abs(landa))
        return "N"
    else
        return "E"
    end
end

end 
