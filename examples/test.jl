using CircuitNetworks, StaticArrays, BenchmarkTools, CUDA, Symbolics

function rlgc_line(r, l, g, c, d)
    function model(f, _)
        z0 = sqrt((r + 2 * π * im * f * l) / (g + 2 * π * im * f * c))
        γ = sqrt((r + 2 * π * im * f * l) * (g + 2 * π * im * f * c))
        s = sinh(γ * d)
        c = cosh(γ * d)
        @SMatrix [c s*z0; s/z0 c]
    end
    #z0 = @SVector [50.0, 50.0]
    #ModelCircuitNetwork(model, z0, Parameter.ABCD)
end

@code_lowered rlgc_line(rand(5)...)

@code_warntype rlgc_line(rand(Float32, 5)...)
