using CircuitNetworks, StaticArrays, BenchmarkTools, CUDA

S = @SMatrix rand(ComplexF32, 2, 2)
z0 = @SArray rand(ComplexF32, 2)
ss = fill(S, 201)
f = range(1e9, 5e9, 201)

n = DataCircuitNetwork(ss, f; kind=Parameter.S)

to_Z(n)