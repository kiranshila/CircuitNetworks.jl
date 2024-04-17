function transducer_gain(S::AbstractMatrix, Γl::Number, Γs::Number)
    ((1 - abs2(Γs)) * abs2(S[2, 1]) * (1 - abs2(Γl))) /
    abs((1 - S[1, 1] * Γs) * (1 - S[2, 2] * Γl) - S[1, 2] * S[2, 1] * Γs * Γl)^2
end

function transducer_gain(n::DataCircuitNetwork{Val{Parameter.S},T,F,Z0}, Γl, Γs) where {T,F,Z0}
    transducer_gain.(n.params, Γl, Γs)
end

"""
    Γout(S,Γs)
Computes the output reflection coefficient of a two port network with its input (Port 1)
terminated with a load whose reflection coefficient is Γs.
"""
function Γout(S::AbstractMatrix, Γs::Number)
    S[2, 2] + (S[1, 2] * S[2, 1] * Γs) / (1 - S[1, 1] * Γs)
end

"""
    Γin(S,Γl)
Computes the input reflection coefficient of a two port network with its output (Port 2)
terminated with a load whose reflection coefficient is Γl.
"""
function Γin(S::AbstractMatrix, Γl::Number)
    S[1, 1] + (S[2, 1] * S[1, 2] * Γl) / (1 - S[2, 2] * Γl)
end

function available_gain(S::AbstractMatrix, Γs::Number)
    # Pozar eqn. 5.88
    pavn = abs2(S[2, 1]) * (1 - abs2(Γs))
    pavs = abs2(1 - S[1, 1] * Γs) * (1 - abs2(Γout(S, Γs)))
    pavn / pavs
end

function available_gain(n::DataCircuitNetwork{Val{Parameter.S},T,F,Z0}, Γs) where {T,F,Z0}
    available_gain.(n.params, Γs)
end

export transducer_gain, available_gain, Γout, Γin