"""The voltage reflection coefficient"""
Γ(ZL::Number, Z0::Number) = (ZL - Z0) / (ZL + Z0)

"""The characteristic impedance of a transmisslion line given the RLGC parameters and frequency"""
Z₀(r::Real, l::Real, g::Real, c::Real, f::Real) = sqrt((r + 2π * im * f * l) / (g + 2π * im * f * c))

"""The propogation constant of a transmisslion line given the RLGC parameters and frequency"""
γ(r::Real, l::Real, g::Real, c::Real, f::Real) = sqrt((r + 2π * im * f * l) * (g + 2π * im * f * c))