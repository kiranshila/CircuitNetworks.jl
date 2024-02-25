using StaticArrays

# We'll use an array of static arrays for the parameters as we can then broadcast over the frequency axis
# This is nice because a Vector of SMatrix is the same layout (just as efficient) as a Matrix but has more
# meaningfull broadcast symantics for CPU

"""
An N-port circuit network whose behavior is defined at a list of associated frequencies

# Fields
- s: The NxNxF S-Parameter matrix
- f: The associated frequency axis for the S-Parameters
- z₀: The characteristic impedance of each port
"""
struct DataCircuitNetwork{S,F,Z0} <: AbstractCircuitNetwork
    s::S
    f::F
    z₀::Z0
    # Write an inner constructor that respects the invariants
    function DataCircuitNetwork{S,F,Z0}(s::S, f::F, z₀::Z0) where {
        T1<:Real,
        T2<:Real,
        T3<:Real,
        SS<:AbstractMatrix{Complex{T1}},
        S<:AbstractVector{SS},
        F<:AbstractVector{T2},
        Z0<:AbstractVector{Complex{T3}},
    }
        for ss in s
            @assert size(ss)[1] == size(ss)[2] "S-Parameters must be square"
            @assert size(ss)[1] == size(first(s))[1] "Every point in frequency must have the same number of ports"
        end
        @assert length(s) == length(f) "Number of S-Parameters must match the number of frequency points"
        @assert size(first(s))[1] == length(z₀) "Number of port impedances must match the number of ports in S-Parameters"
        new(s, f, z₀)
    end
end

# Manual outer constructor
function DataCircuitNetwork(s::S, f::F, z₀::Z0) where {
    T1<:Real,
    T2<:Real,
    T3<:Real,
    SS<:AbstractMatrix{Complex{T1}},
    S<:AbstractVector{SS},
    F<:AbstractVector{T2},
    Z0<:AbstractVector{Complex{T3}},
}
    DataCircuitNetwork{S,F,Z0}(s, f, z₀)
end

function show(io::IO, network::DataCircuitNetwork)
    N = length(network.z₀)
    f_label, f_scale = freq_scale(network)
    print(io, "$(N)-Port Network: ")
    print(io, "$(minimum(network.f)/f_scale)-$(maximum(network.f)/f_scale) $(f_label), ")
    print(io, "$(length(network.f)) pts, ")
    print(io, "z₀=[$(join(network.z₀,", "))]")
end

# Convenience constructors

# Real vector of port impedance
function DataCircuitNetwork(s::S, f::F, z₀::Z0) where {
    T1<:Real,
    T2<:Real,
    T3<:Real,
    SS<:AbstractMatrix{Complex{T1}},
    S<:AbstractVector{SS},
    F<:AbstractVector{T2},
    Z0<:AbstractVector{T3},
}
    DataCircuitNetwork(s, f, convert.(Complex, z₀))
end

# Same port impedance
function DataCircuitNetwork(s::AbstractVector{SS}, f::AbstractVector{T2}, z₀::Number) where {
    T1<:Real,
    T2<:Real,
    SS<:AbstractMatrix{Complex{T1}},
}
    DataCircuitNetwork(s, f, fill(z₀, size(first(s))[1]))
end

# Default port impedance (50Ω)
function DataCircuitNetwork(s::AbstractVector{SS}, f::AbstractVector{T2}) where {
    T1<:Real,
    T2<:Real,
    SS<:AbstractMatrix{Complex{T1}},
}
    DataCircuitNetwork(s, f, fill(50.0, size(first(s))[1]))
end

# Public exports
export DataCircuitNetwork