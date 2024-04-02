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

    @testset "Coax" begin
        # Double checking with linecalc results
        # Not sure the method of linecalc, so the results may be slightly different

        freq = 10e9
        d_i = 0.5e-3
        d_o = 1.5e-3
        d = 3e-3

        # Copper/Teflon Coax
        σc = 5.813e7
        εr = 2.1
        tanδ = 0.002

        r, l, g, c = coax_rlgc(d_i / 2, d_o / 2, freq; εr=εr, tanδ=tanδ, σc=σc)

        z0 = rlgc_z0(r, l, g, c, freq)
        γ = rlgc_γ(r, l, g, c, freq)

        @test real(z0) ≈ 45.4554 atol = 0.01

        # α -> neppers/m to db/m
        a_db = 0.014
        α = a_db * 0.1151277918 / d
        @test real(γ) ≈ α atol = 0.1

        # β is 2π/λ
        λ = (299792458 / freq) * (1 / √(εr))
        β = 2π / λ
        @test imag(γ) ≈ β atol = 0.01

        abcd = abcd_tline(γ, z0, d)
        s = a2s(abcd)
        @test angle(s[2, 1]) * 180 / π ≈ -52.2051 atol = 0.2 # Why is this one so off?
    end

    @testset "Noise Correctness" begin
        freq = 1e9
        rlgc = coax_rlgc(0.9e-3 / 2, 3.275e-3 / 2, freq; εr=2.1, tanδ=0.002, σc=5.813e7)
        abcd = rlgc2abcd(rlgc..., 25e-3, freq)
        s = a2s(abcd) # 50 Ohm default
        nf = noise_figure(s, 0, 290) # Matched generator to 50 Ohms, 290K physical temperature
        nt = noise_temperature(nf)
        # Compare to ADS
        @test nt ≈ 0.953 atol = 0.01
    end
end
