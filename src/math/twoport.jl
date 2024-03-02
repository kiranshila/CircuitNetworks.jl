using Tullio

function transducer_gain(S::AbstractMatrix, Γl::Number, Γs::Number)
    ((1 - abs(Γs)^2) * abs(S[2, 1])^2 * (1 - abs(Γl)^2)) /
    abs((1 - S[1, 1] * Γs) * (1 - S[2, 2] * Γl) - S[1, 2] * S[2, 1] * Γs * Γl)^2
end

function transducer_gain(s::T, Γl::AbstractVector, Γs::AbstractVector) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    @tullio grad = Dual g[i] := transducer_gain(s[i], Γl[i], Γs[i])
end

function transducer_gain(s::T, Γl::Number, Γs::AbstractVector) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    @tullio grad = Dual g[i] := transducer_gain(s[i], Γl, Γs[i])
end

function transducer_gain(s::T, Γl::AbstractVector, Γs::Number) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    @tullio grad = Dual g[i] := transducer_gain(s[i], Γl[i], Γs)
end

function transducer_gain(s::T, Γl::Number, Γs::Number) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    @tullio grad = Dual g[i] := transducer_gain(s[i], Γl, Γs)
end

function transducer_gain(n::DataCircuitNetwork{Val{Parameter.S},T,F,Z0}, Γl, Γs) where {T,F,Z0}
    transducer_gain(n.params, Γl, Γs)
end

function available_gain(S::AbstractMatrix, Γs::Number)
    Γout = S[2, 2] + (S[1, 2] * S[2, 1] * Γs) / (1 - S[1, 1] * Γs)
    ((1 - abs(Γs)^2) * abs(S[2, 1])^2) /
    (abs(1 - S[1, 1] * Γs)^2 * (1 - abs(Γout)^2))
end

function available_gain(s::T, Γs::AbstractVector) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    @tullio grad = Dual g[i] := available_gain(s[i], Γs[i])
end

function available_gain(s::T, Γs::Number) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    @tullio grad = Dual g[i] := available_gain(s[i], Γs)
end

function available_gain(n::DataCircuitNetwork{Val{Parameter.S},T,F,Z0}, Γs) where {T,F,Z0}
    available_gain(n.params, Γs)
end

export transducer_gain, available_gain