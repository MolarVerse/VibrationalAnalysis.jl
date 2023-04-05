using VibrationalAnalysis
using Test

# Test symmetrize matrices src/symmetrize.jl

@testset "Symmetrize Matrices" begin

    # Test Symmetrize Addition
    @test symmetrize_addition([2.0 3.0; 1.0 -2.0]) == [2.0 2.0; 2.0 -2.0]

    # Test Symmetrize Multiplication
    @test symmetrize_multiplication([2.0 3.0; 1.0 -2.0]) == [2.1213203435596424 0.7071067811865476; 0.7071067811865475 3.5355339059327373]

end

@testset "Mass-weighted hessian" begin
        # Test Mass-weighted hessian

        hessian = read_hessian("../data/test_hessian_h2o.dat")
        _, atom_masses , _ = read_rst("../data/test_h2o.rst")
        
        mass_weight_hessian = mass_weighted_hessian(hessian, atom_masses)
        @test mass_weight_hessian isa Matrix{Float64}

        # Test Symmetrize Multiplication
        hessian = [4.0 2.0 0.5; 2.0 4.0 0.0; 0.5 0.0 4.0]
        atom_masses = [2.0]
        @test mass_weighted_hessian(hessian, atom_masses) ≈ [2.0 1.0 0.25; 1.0 2.0 0.0; 0.25 0.0 2.0] atol=1e-5
        @test mass_weighted_hessian_add(hessian, atom_masses) == [2.0 1.0 0.25; 1.0 2.0 0.0; 0.25 0.0 2.0] 

end