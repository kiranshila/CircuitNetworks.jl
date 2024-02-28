using CircuitNetworks, Test, StaticArrays

@testset "CircuitNetworks.jl" begin

    @testset "Conversions" begin
        f_n = 10
        for n in [1, 2, 3, 5]
            z0 = @SVector rand(ComplexF64, n)
            m = [@SMatrix rand(ComplexF64, n, n) for _ in 1:f_n]
            @test z2s(s2z(m, z0), z0) ≈ m
            @test y2s(s2y(m, z0), z0) ≈ m
            @test z2y(y2z(m)) ≈ m

            # Test 2-port only networks
            if n == 2
                @test a2s(s2a(m, z0), z0) ≈ m
            end
        end
    end

    @testset "Power Gains" begin
        m = @SMatrix [
            +0.8510+0.1493im -0.1250+0.4847im;
            -0.1250+0.4847im +0.8180+0.2771im
        ]
        z0 = 50 + 0im
        zs = 50 + 0im
        zl = 50 + 0im

    end

end
