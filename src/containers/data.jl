using StaticArrays
import Base.show

# We'll use an array of static arrays for the parameters as we can then broadcast over the frequency axis
# This is nice because a Vector of SMatrix is the same layout (just as efficient) as a Matrix but has more
# meaningfull broadcast symantics for CPU

"""
An N-port circuit network whose behavior is defined at a list of associated frequencies

# Fields
- params: The NxNxF port parameter matrix
- f: The associated frequency axis for the parameters
- z₀: The characteristic impedance of each port
- kind: The representation of the port parameters
"""
struct DataCircuitNetwork{P,T,F,Z0} <: AbstractCircuitNetwork
    params::T
    f::F
    z₀::Z0
    # Write an inner constructor that respects the invariants
    function DataCircuitNetwork{P,T,F,Z0}(params::T, f::F, z₀::Z0) where {
        T1<:Real,
        T2<:Real,
        T3<:Real,
        TT<:AbstractMatrix{Complex{T1}},
        T<:AbstractVector{TT},
        F<:AbstractVector{T2},
        Z0<:AbstractVector{Complex{T3}},
        P
    }
        for p in params
            @assert size(p)[1] == size(p)[2] "Parameters must be square"
            @assert size(p)[1] == size(first(params))[1] "Every point in frequency must have the same number of ports"
        end
        @assert length(params) == length(f) "Number of parameters must match the number of frequency points"
        @assert size(first(params))[1] == length(z₀) "Number of port impedances must match the number of ports in the parameters"
        new(params, f, z₀)
    end
end

# Manual outer constructor
function DataCircuitNetwork(params::T, f::F, z₀::Z0, kind::ParameterType) where {
    T1<:Real,
    T2<:Real,
    T3<:Real,
    TT<:AbstractMatrix{Complex{T1}},
    T<:AbstractVector{TT},
    F<:AbstractVector{T2},
    Z0<:AbstractVector{Complex{T3}},
}
    DataCircuitNetwork{Val{kind},T,F,Z0}(params, f, z₀)
end

function show(io::IO, network::DataCircuitNetwork{Val{P},T,F,Z0}) where {P,T,F,Z0}
    N = length(network.z₀)
    f_label, f_scale = freq_scale(network)
    t = String(P)
    print(io, "$(N)-Port $(t)-Parameter Network: ")
    print(io, "$(minimum(network.f)/f_scale)-$(maximum(network.f)/f_scale) $(f_label), ")
    print(io, "$(length(network.f)) pts, ")
    print(io, "z₀=[$(join(network.z₀,", "))]")
end

# Convenience constructors

# Real vector of port impedance
function DataCircuitNetwork(params::T, f::F, z₀::Z0, kind::ParameterType) where {
    T1<:Real,
    T2<:Real,
    T3<:Real,
    TT<:AbstractMatrix{Complex{T1}},
    T<:AbstractVector{TT},
    F<:AbstractVector{T2},
    Z0<:AbstractVector{T3},
}
    DataCircuitNetwork(params, f, convert.(Complex, z₀), kind)
end

# Same port impedance for both ports, as kwargs with defaults
function DataCircuitNetwork(params::AbstractVector{TT}, f::AbstractVector{T2}; z₀::Number=50.0, kind::ParameterType=Parameter.S) where {
    T1<:Real,
    T2<:Real,
    TT<:AbstractMatrix{Complex{T1}},
}
    DataCircuitNetwork(params, f, fill(z₀, size(first(params))[1]), kind)
end

# Public exports
export DataCircuitNetwork