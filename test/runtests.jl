using CircuitNetworks, Test, StaticArrays

@testset "CircuitNetworks.jl" begin

    @testset "Conversion Round Trip" begin
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

    @testset "Conversion Correctness" begin
        # Following MATLAB (https://www.mathworks.com/help/rf/ref/s2z.html)
        s_11 = 0.61 * cispi(165 / 180)
        s_21 = 3.72 * cispi(59 / 180)
        s_12 = 0.05 * cispi(42 / 180)
        s_22 = 0.45 * cispi(-48 / 180)
        s = [s_11 s_12; s_21 s_22]
        @testset "S<->Z" begin
            @test s2z(s) ≈ [0.1141+0.1567im 0.0352+0.0209im; 2.0461+2.2524im 0.7498-0.3803im] .* 1e2 atol = 0.01
        end
        @testset "S<->Y" begin
            @test s2y(s) ≈ [0.0647-0.0059im -0.0019-0.0025im; -0.0826-0.2200im 0.0037+0.0145im] atol = 0.01
        end
        @testset "S<->ABCD" begin
            @test s2a(s) ≈ [0.0633+0.0069im 1.4958-3.9839im; 0.0022-0.0024im 0.0732-0.2664im] atol = 0.01
        end
    end

    @testset "Power Gain Correctness" begin
        # Following MATLAB (https://www.mathworks.com/help/rf/ref/powergain.html)
        s_11 = 0.61 * cispi(165 / 180)
        s_21 = 3.72 * cispi(59 / 180)
        s_12 = 0.05 * cispi(42 / 180)
        s_22 = 0.45 * cispi(-48 / 180)
        s = [s_11 s_12; s_21 s_22]
        z0 = 50
        zs = 10 + 20im
        zl = 30 - 40im

        Γl = CircuitNetworks.Γ(zl, z0)
        Γs = CircuitNetworks.Γ(zs, z0)

        @test transducer_gain(s, Γl, Γs) ≈ 4.7066 atol = 0.01
        @test available_gain(s, Γs) ≈ 11.4361 atol = 0.01
    end

    @testset "Noise Correctness" begin end
end
