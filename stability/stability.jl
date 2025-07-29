module Stability    

export stability

using LinearAlgebra
using IntervalArithmetic
using RadiiPolynomial

function stability(landa)

    if isstrictless(interval(1), abs(landa))
        return "U"
    elseif in_interval(1, abs(landa))
        return "N"
    else
        return "S"
    end
end

end 
