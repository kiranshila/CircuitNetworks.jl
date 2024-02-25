"""
Utility function to get the frequency scale across a range for printing and plotting
"""
function freq_scale(freqs::AbstractArray{T}) where {T<:Real}
    max_f = maximum(freqs)
    if max_f < 1e3
        ("Hz", 1)
    elseif max_f < 1e6
        ("kHz", 1e3)
    elseif max_f < 1e9
        ("MHz", 1e6)
    elseif max_f < 1e12
        ("GHz", 1e9)
    elseif max_f < 1e15
        ("THz", 1e12)
    else
        ("Hz", 1)
    end
end

freq_scale(network::DataCircuitNetwork) = freq_scale(network.f)