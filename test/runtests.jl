using AnalyticalEOR
using Test

@testset "AnalyticalEOR.jl" begin

    # Test relative permeability functions
    @test water_rel_perm(0.5, 0.2, 0.2, 0.3, 3.0) ≈ 0.0374999999
    @test water_rel_perm([0.5, 0.5], 0.2, 0.2, 0.3, 3.0) ≈ Vector([0.0375, 0.0375])
   
    @test oil_rel_perm(0.5, 0.2, 0.2, 0.3, 3.0) ≈ 0.0375
    @test oil_rel_perm([0.5, 0.5], 0.2, 0.2, 0.3, 3.0) ≈ Vector([0.0375, 0.0375])

    @test fractional_flow([0.5, 0.6], 0.2, 0.2, 0.3, 1.0, 3.0, 3.0, 1.0, 1.0) ≈ Vector([ 0.2307692307692305,  0.7058823529411761])
end

@testset "Derivative Tests" begin
    
    @test krw_derivative([0.5, 0.5], 0.2, 0.2, 0.3, 3.0) ≈ [0.375, 0.375]
    @test kro_derivative([0.5, 0.5], 0.2, 0.2, 0.3, 3.0) ≈ [-0.375, -0.375]

    # test derivative of fractional flow
end

