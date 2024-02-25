dB(x::Real) = 10 * log10(x)
dB(x::Complex) = 10 * log10(abs(x))

"""The voltage reflection coefficient"""
Γ(Z::Number, Z0::Number) = (Z - Z0) / (Z + Z0)