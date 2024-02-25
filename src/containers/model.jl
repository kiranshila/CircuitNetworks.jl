"""
An N-port circuit network whose behavior is defined functionally.

# Fields
 - f: A function that accepts a frequency and N-long vector of port impedances and returns an NxN S-Matrix
 - z₀: A N-long vector of port impedances for the modeled NxN network
"""
struct ModeledCircuitNetwork{Z0,F} <: AbstractCircuitNetwork
    f::F
    z₀::Z0
end

# TODO Create constructor for ModeledCircuitNetwork that checks invariants