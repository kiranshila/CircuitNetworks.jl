function noise_figure(nfmin, rn, Γopt, Γs; z0=50.0)
    nfmin + (4 * rn * (abs(Γs - Γopt)^2)) / (z0 * (1 - abs(Γs)^2) * abs(1 + Γopt)^2)
end

export noise_figure