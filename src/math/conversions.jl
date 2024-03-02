using LinearAlgebra, Tullio

#### S <-> Z

# Power Waves and the Scattering Network (paper)

function s2z(S::AbstractMatrix, z0::AbstractVector)
    F = Diagonal(@. 1 / (2 * sqrt(abs(real(z0)))))
    G = Diagonal(z0)
    inv(F) * inv(I - S) * (S * G + G') * F
end

function s2z(S::AbstractMatrix; z0::Number=50.0)
    s2z(S,fill(z0,size(S)[1]))
end

function s2z!(z::T, s::T, z0::AbstractVector) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    @tullio grad=Dual z[i] = s2z(s[i], z0)
end

function s2z(s::T, z0::AbstractVector) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    z = similar(s)
    s2z!(z, s, z0)
end

function z2s(Z::AbstractMatrix, z0::AbstractVector)
    F = Diagonal(@. 1 / (2 * sqrt(abs(real(z0)))))
    G = Diagonal(z0)
    F * (Z - G') * inv(Z + G) * inv(F)
end

function z2s(Z::AbstractMatrix; z0::Number=50.0)
    z2s(Z,fill(z0,size(Z)[1]))
end

function z2s!(s::T, z::T, z0::AbstractVector) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    @tullio grad=Dual s[i] = z2s(z[i], z0)
end

function z2s(z::T, z0::AbstractVector) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    s = similar(z)
    z2s!(s, z, z0)
end

#### S <-> Y

function s2y(S::AbstractMatrix, z0::AbstractVector)
    F = Diagonal(@. 1 / (2 * sqrt(abs(real(z0)))))
    G = Diagonal(z0)
    inv(F) * inv(S * G + G') * (I - S) * F
end

function s2y(S::AbstractMatrix; z0::Number=50.0)
    s2y(S,fill(z0,size(S)[1]))
end

function s2y!(y::T, s::T, z0::AbstractVector) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    @tullio grad=Dual y[i] = s2y(s[i], z0)
end

function s2y(s::T, z0::AbstractVector) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    y = similar(s)
    s2y!(y, s, z0)
end

function y2s(Y::AbstractMatrix, z0::AbstractVector)
    F = Diagonal(@. 1 / (2 * sqrt(abs(real(z0)))))
    G = Diagonal(z0)
    F * (I - G' * Y) * inv(I + G * Y) * inv(F)
end

function y2s(Y::AbstractMatrix; z0::Number=50.0)
    y2s(Y,fill(z0,size(S)[1]))
end

function y2s!(s::T, y::T, z0::AbstractVector) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    @tullio grad=Dual s[i] = y2s(y[i], z0)
end

function y2s(y::T, z0::AbstractVector) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    s = similar(y)
    y2s!(s, y, z0)
end

#### Y <-> Z

function y2z!(z::T, y::T) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    @tullio grad=Dual z[i] = inv(y[i])
end

function y2z(y::T) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    z = similar(y)
    y2z!(z, y)
end

function z2y!(y::T, z::T) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    @tullio grad=Dual y[i] = inv(z[i])
end

function z2y(z::T) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    y = similar(z)
    z2y!(y, z)
end

#### S <-> ABCD (2-Port only)

# Conversions between s,z,y,h, abcd, and t which are valid for complex source and load impedances

function s2a(S::AbstractMatrix, z0::AbstractVector)
    denom = 2 * S[2, 1] * sqrt(real(z0[1]) * real(z0[2]))
    @SMatrix [
        ((z0[1]' + S[1, 1] * z0[1]) * (1 - S[2, 2]) + S[1, 2] * S[2, 1] * z0[1]) / denom#=
        =# ((z0[1]' + S[1, 1] * z0[1]) * (z0[2]' + S[2, 2] * z0[2]) - S[1, 2] * S[2, 1] * z0[1] * z0[2]) / denom;
        ((1 - S[1, 1]) * (1 - S[2, 2]) - S[1, 2] * S[2, 1]) / denom#=
        =# ((1 - S[1, 1]) * (z0[2]' + S[2, 2] * z0[2]) + S[1, 2] * S[2, 1] * z0[2]) / denom
    ]
end

function s2a(S::AbstractMatrix; z0::Number=50.0)
    s2a(S,fill(z0,size(S)[1]))
end

function s2a!(a::T, s::T, z0::AbstractVector) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    @tullio grad=Dual a[i] = s2a(s[i], z0)
end

function s2a(s::T, z0::AbstractVector) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    a = similar(s)
    s2a!(a, s, z0)
end

function a2s(A::AbstractMatrix, z0::AbstractVector)
    denom = A[1, 1] * z0[2] + A[1, 2] + A[2, 1] * z0[1] * z0[2] + A[2, 2] * z0[1]
    @SMatrix [
        (A[1, 1] * z0[2] + A[1, 2] - A[2, 1] * z0[1]' * z0[2] - A[2, 2] * z0[1]') / denom#=
        =#(2 * (A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1]) * sqrt(real(z0[1] * real(z0[2])))) / denom;
        (2 * sqrt(real(z0[1]) * real(z0[2]))) / denom#=
        =#(-A[1, 1] * z0[2]' + A[1, 2] - A[2, 1] * z0[1] * z0[2]' + A[2, 2] * z0[1]) / denom
    ]
end

function a2s(A::AbstractMatrix; z0::Number=50.0)
    a2s(A,fill(z0,size(A)[1]))
end

function a2s!(s::T, a::T, z0::AbstractVector) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    @tullio grad=Dual s[i] = a2s(a[i], z0)
end

function a2s(a::T, z0::AbstractVector) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    s = similar(a)
    a2s!(s, a, z0)
end

# Public exports
export
    s2z!, s2z, z2s!, z2s,
    s2y!, s2y, y2s!, y2s,
    y2z!, y2z, z2y!, z2y,
    s2a!, s2a, a2s!, a2s

#############################  Network Conversions ################################

# Multiple dispatch is handy here

"""
Construct a new `DataCircuitNetwork` by converting another to S parameters
"""
to_S(n::DataCircuitNetwork{Val{Parameter.S}, T, F, Z0}) where {T,F,Z0} = n

function to_S(n::DataCircuitNetwork{Val{Parameter.Z}, T, F, Z0}) where {T,F,Z0}
    params = z2s(n.params, n.z₀)
    DataCircuitNetwork(params,n.f, n.z₀, Parameter.S)
end

function to_S(n::DataCircuitNetwork{Val{Parameter.Y}, T, F, Z0}) where {T,F,Z0}
    params = y2s(n.params, n.z₀)
    DataCircuitNetwork(params,n.f, n.z₀, Parameter.S)
end

function to_S(n::DataCircuitNetwork{Val{Parameter.ABCD}, T, F, Z0}) where {T,F,Z0}
    params = a2s(n.params, n.z₀)
    DataCircuitNetwork(params,n.f, n.z₀, Parameter.S)
end

"""
Construct a new `DataCircuitNetwork` by converting another to ABCD parameters
"""
to_ABCD(n::DataCircuitNetwork{Val{Parameter.ABCD}, T, F, Z0}) where {T,F,Z0} = n

function to_ABCD(n::DataCircuitNetwork{Val{Parameter.S}, T, F, Z0}) where {T,F,Z0}
    params = s2a(n.params, n.z₀)
    DataCircuitNetwork(params,n.f, n.z₀, Parameter.ABCD)
end

"""
Construct a new `DataCircuitNetwork` by converting another to Z parameters
"""
to_Z(n::DataCircuitNetwork{Val{Parameter.Z}, T, F, Z0}) where {T,F,Z0} = n

function to_Z(n::DataCircuitNetwork{Val{Parameter.S}, T, F, Z0}) where {T,F,Z0}
    params = s2z(n.params, n.z₀)
    DataCircuitNetwork(params,n.f, n.z₀, Parameter.Z)
end

"""
Construct a new `DataCircuitNetwork` by converting another to Y parameters
"""
to_Y(n::DataCircuitNetwork{Val{Parameter.Y}, T, F, Z0}) where {T,F,Z0} = n

function to_Y(n::DataCircuitNetwork{Val{Parameter.S}, T, F, Z0}) where {T,F,Z0}
    params = s2y(n.params, n.z₀)
    DataCircuitNetwork(params,n.f, n.z₀, Parameter.Y)
end

export to_S, to_ABCD, to_Y, to_Z