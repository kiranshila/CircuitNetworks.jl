const T₀ = 290 # K, IEEE Standard Reference Temperature

"""
Noise figure from the four noise paramters and a source reflection Γs (linear)
"""
function noise_figure(nfmin::Number, rn::Number, Γopt::Number, Γs::Number; z0=50.0)
    nfmin + (4 * rn * (abs2(Γs - Γopt))) / (z0 * (1 - abs2(Γs)) * abs2(1 + Γopt))
end

"""
Noise figure of a passive two port network (linear)
"""
function noise_figure(S::AbstractMatrix, Γs::Number, T::Number)
    ga = available_gain(S, Γs)
    1 + ((1 - ga) / ga) * (T / T₀)
end

"""
Convert noise figure (linear) to noise temperature
"""
noise_temperature(nf::Number) = T₀ * (nf - 1)

export noise_figure, noise_temperature