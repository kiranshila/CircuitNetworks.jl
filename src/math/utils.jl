dB(x::Real) = 10 * log10(x)
dB(x::Complex) = 10 * log10(abs(x))

"""The voltage reflection coefficient"""
Γ(Z::Number, Z0::Number) = (Z - Z0) / (Z + Z0)

"""The characteristic impedance of a transmisslion line given the RLGC parameters and frequency"""
Z₀(r::Real, l::Real, g::Real, c::Real, f::Real) = sqrt((r + 2π * im * f * l) / (g + 2π * im * f * c))

"""The propogation constant of a transmisslion line given the RLGC parameters and frequency"""
γ(r::Real, l::Real, g::Real, c::Real, f::Real) = sqrt((r + 2π * im * f * l) * (g + 2π * im * f * c))