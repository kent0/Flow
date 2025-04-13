using Metal, LinearAlgebra, BenchmarkTools

#n = p = m = 512

n1d=512

n = n1d * n1d
p = n1d
m = n1d

flops = n*m*(2p-1)

a = MtlArray(rand(Float32, n, p));
b = MtlArray(rand(Float32, p, m));
c = MtlArray(zeros(Float32, n, m));

bench = @benchmark Metal.@sync mul!(c, a, b)
gflops = flops / (mean(bench.times)/1e9)
println("Metal GFLOPS: ", gflops)

# Initialize matrices with random values and zeros
a = rand(Float32, n, p)
b = rand(Float32, p, m)
c = zeros(Float32, n, m)

# Benchmark the matrix multiplication
bench = @benchmark mul!(c, a, b)
gflops = flops / (mean(bench.times) / 1e9)

println("GFLOPS: ", gflops)
