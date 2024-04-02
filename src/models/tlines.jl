# Transmission line models and equations (TEM)

const ε₀ = 8.85418782e-12
const μ₀ = 1.25663706e-6

"""
Compute the RLGC (per unit length transmission line parameters) of coax

## Geometric Parameters
- `a`: Radius of the inner conductor
- `b`: Radius of the outer conductor
## Electrical Parameters
- `freq`: Frequency of operation
## Conductor Parameters
- `σc`: Conductivity of the conductors (S/m)
- `μr_c`:  Relative permeability of the conductors
## Dielectric Parameters
- `μr_d`:  Relative permeability of the insulating material between conductors
- `εr`: Relative permitivity of the insulating material between conductors
- `tanδ`: Loss tangent of the dielectric
"""
function coax_rlgc(a, b, freq; μr_d=1, εr=1, μr_c=1, tanδ=0, σc=Inf)
    μc = μr_c * μ₀
    μd = μr_d * μ₀
    εd = εr * ε₀

    ε″ = εd * tanδ
    lba = log(b / a)
    ω = 2π * freq

    R = if isinf(σc)
        0.0
    else
        sqrt(π * freq * μc * σc) / (2π * σc) * (1 / a + 1 / b)
    end
    L = μd / (2π) * lba
    G = 2π * ω * ε″ / lba
    C = (2π * εd) / lba

    (R, L, G, C)
end

"""The characteristic impedance of a transmission line given the RLGC parameters and frequency"""
function rlgc_z0(r::Real, l::Real, g::Real, c::Real, freq::Real)
    sqrt((r + 2π * im * freq * l) / (g + 2π * im * freq * c))
end

"""The propogation constant of a transmission line given the RLGC parameters and frequency"""
function rlgc_γ(r::Real, l::Real, g::Real, c::Real, freq::Real)
    sqrt((r + 2π * im * freq * l) * (g + 2π * im * freq * c))
end

"""
Compute the ABCD matrix of a trasmission line by its complex propogation constant `γ`, its characteristic impedance `z₀` and the length `d`.
"""
function abcd_tline(γ, z₀, d)
    sh = sinh(γ * d)
    ch = cosh(γ * d)
    @SMatrix [ch sh*z₀; sh/z₀ ch]
end

"""
Construct an ABCD matrix of a transmission line from its RLGC paratemeters, the length of the line `d` and the frequency of operation `freq`
"""
rlgc2abcd(r, l, g, c, d, freq) = abcd_tline(rlgc_γ(r, l, g, c, freq), rlgc_z0(r, l, g, c, freq), d)

export coax_rlgc, rlgc2abcd, abcd_tline, rlgc2abcd, rlgc_z0, rlgc_γ