# Indexing operations for an N-port network

import Base.getindex, Base.view

# Parameter access

function getindex(n::DataCircuitNetwork, kind::AbstractString)
    # Conversions between network types
    if lowercase(kind) == "s"
        n.s
    elseif lowercase(kind) == "abcd"
        s2a(n.s, n.z₀)
    elseif lowercase(kind) == "z"
        s2z(n.s, n.z₀)
    elseif lowercase(kind) == "y"
        s2y(n.s, n.z₀)
    else
        error("$kind is not a valid parameter type")
    end
end

# We'll have to be careful here being correct and avoiding copies when indexing down the frequency axis