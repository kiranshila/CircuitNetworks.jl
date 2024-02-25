# There might be a better way to formulate this that interacts with ModellingToolkit better

"""
An N-port circuit network whose behavior is defined functionally.

# Fields
 - f: A function that accepts a frequency and N-long vector of port impedances and returns an NxN port parameter matrix
 - z₀: A N-long vector of port impedances for the modeled NxN network
"""
struct ModelCircuitNetwork{P,Z0,F} <: AbstractCircuitNetwork
    f::F
    z₀::Z0
    function ModelCircuitNetwork{P,Z0,F}(f::F, z0::Z0) where {P,F,Z0}
        new(f, z0)
    end
end

function ModelCircuitNetwork(f::F, z0::Z0, kind::ParameterType) where {F<:Function,Z0<:AbstractVector}
    ModelCircuitNetwork{Val{kind},Z0,F}(f, z0)
end

(n::ModelCircuitNetwork)(f) = n.f(f, n.z₀)

# TODO Create constructor for ModeledCircuitNetwork that checks invariants

export ModelCircuitNetwork