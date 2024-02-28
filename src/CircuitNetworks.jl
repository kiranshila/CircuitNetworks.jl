module CircuitNetworks

using EnumX

"""The abstract type for all circuit networks"""
abstract type AbstractCircuitNetwork end

"""
Paramater types representable in `CircuitNetworks`
"""
@enumx T = ParameterType Parameter begin
    "S-parameters"
    S
    "Z-parameters"
    Z
    "Y-parameters"
    Y
    "ABCD-parameters"
    ABCD
end
import .Parameter.ParameterType
export Parameter, ParameterType

function ParameterType(p::AbstractString)
    p = lowercase(p)
    if p == "s"
        Parameter.S
    elseif p == "z"
        Parameter.Z
    elseif p == "Y"
        Parameter.Y
    elseif p == "ABCD"
        Parameter.ABCD
    else
        error("No matching parameter type")
    end
end

include("containers/data.jl")
include("containers/model.jl")
include("visualizations/utils.jl")
include("math/utils.jl")
include("math/conversions.jl")
include("math/cascade.jl")
include("math/twoport.jl")
include("math/noise.jl")
include("containers/introspection.jl")

end
