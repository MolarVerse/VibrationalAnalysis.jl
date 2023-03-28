using VibrationalAnalysis
using Test

# Test symmetrize matrices src/symmetrize.jl

@testset "Symmetrize Matrices" begin

    # Test Symmetrize Addition
    @test symmetrize_addition([2.0 3.0; 1.0 -2.0]) == [2.0 2.0; 2.0 -2.0]

    # Test Symmetrize Multiplication
    @test symmetrize_multiplication([2.0 3.0; 1.0 -2.0]) == [2.1213203435596424 0.7071067811865476; 0.7071067811865475 3.5355339059327373]

end