using LinearAlgebra, Tullio

#### S <-> Z

# Power Waves and the Scattering Network (paper)

function s2z_kern(S::AbstractMatrix, z0::AbstractVector)
    F = Diagonal(@. 1 / (2 * sqrt(abs(real(z0)))))
    G = Diagonal(z0)
    inv(F) * inv(I - S) * (S * G + G') * F
end

function s2z!(z::T, s::T, z0::AbstractVector) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    @tullio z[i] = s2z_kern(s[i], z0)
end

function s2z(s::T, z0::AbstractVector) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    z = similar(s)
    s2z!(z, s, z0)
end

function z2s_kern(Z::AbstractMatrix, z0::AbstractVector)
    F = Diagonal(@. 1 / (2 * sqrt(abs(real(z0)))))
    G = Diagonal(z0)
    F * (Z - G') * inv(Z + G) * inv(F)
end

function z2s!(s::T, z::T, z0::AbstractVector) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    @tullio s[i] = z2s_kern(z[i], z0)
end

function z2s(z::T, z0::AbstractVector) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    s = similar(z)
    z2s!(s, z, z0)
end

#### S <-> Y

function s2y_kern(S::AbstractMatrix, z0::AbstractVector)
    F = Diagonal(@. 1 / (2 * sqrt(abs(real(z0)))))
    G = Diagonal(z0)
    inv(F) * inv(S * G + G') * (I - S) * F
end

function s2y!(y::T, s::T, z0::AbstractVector) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    @tullio y[i] = s2y_kern(s[i], z0)
end

function s2y(s::T, z0::AbstractVector) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    y = similar(s)
    s2y!(y, s, z0)
end

function y2s_kern(Y::AbstractMatrix, z0::AbstractVector)
    F = Diagonal(@. 1 / (2 * sqrt(abs(real(z0)))))
    G = Diagonal(z0)
    F * (I - G' * Y) * inv(I + G * Y) * inv(F)
end

function y2s!(s::T, y::T, z0::AbstractVector) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    @tullio s[i] = y2s_kern(y[i], z0)
end

function y2s(y::T, z0::AbstractVector) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    s = similar(y)
    y2s!(s, y, z0)
end

#### Y <-> Z

function y2z!(z::T, y::T) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    @tullio z[i] = inv(y[i])
end

function y2z(y::T) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    z = similar(y)
    y2z!(z, y)
end

function z2y!(y::T, z::T) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    @tullio y[i] = inv(z[i])
end

function z2y(z::T) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    y = similar(z)
    z2y!(y, z)
end

#### S <-> ABCD (2-Port only)

# Conversions between s,z,y,h, abcd, and t which are valid for complex source and load impedances

function s2a_kern(S::AbstractMatrix, z0::AbstractVector)
    denom = 2 * S[2, 1] * sqrt(real(z0[1]) * real(z0[2]))
    @SMatrix [
        ((z0[1]' + S[1, 1] * z0[1]) * (1 - S[2, 2]) + S[1, 2] * S[2, 1] * z0[1]) / denom#=
        =# ((z0[1]' + S[1, 1] * z0[1]) * (z0[2]' + S[2, 2] * z0[2]) - S[1, 2] * S[2, 1] * z0[1] * z0[2]) / denom;
        ((1 - S[1, 1]) * (1 - S[2, 2]) - S[1, 2] * S[2, 1]) / denom#=
        =# ((1 - S[1, 1]) * (z0[2]' + S[2, 2] * z0[2]) + S[1, 2] * S[2, 1] * z0[2]) / denom
    ]
end

function s2a!(a::T, s::T, z0::AbstractVector) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    @tullio a[i] = s2a_kern(s[i], z0)
end

function s2a(s::T, z0::AbstractVector) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    a = similar(s)
    s2a!(a, s, z0)
end

function a2s_kern(A::AbstractMatrix, z0::AbstractVector)
    denom = A[1, 1] * z0[2] + A[1, 2] + A[2, 1] * z0[1] * z0[2] + A[2, 2] * z0[1]
    @SMatrix [
        (A[1, 1] * z0[2] + A[1, 2] - A[2, 1] * z0[1]' * z0[2] - A[2, 2] * z0[1]') / denom#=
        =#(2 * (A[1, 1] * A[2, 2] - A[1, 2] * A[2, 1]) * sqrt(real(z0[1] * real(z0[2])))) / denom;
        (2 * sqrt(real(z0[1]) * real(z0[2]))) / denom#=
        =#(-A[1, 1] * z0[2]' + A[1, 2] - A[2, 1] * z0[1] * z0[2]' + A[2, 2] * z0[1]) / denom
    ]
end

function a2s!(s::T, a::T, z0::AbstractVector) where {TT<:AbstractMatrix,T<:AbstractVector{TT}}
    @tullio s[i] = a2s_kern(a[i], z0)
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