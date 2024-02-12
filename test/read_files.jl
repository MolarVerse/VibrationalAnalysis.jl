using VibrationalAnalysis
using Test

# Test reading RST files src/read_files.jl

@testset "Read Restart File Exceptions" begin
    @test_throws ErrorException read_rst("data/not_a_file.rst")
    @test_throws ErrorException read_rst("data/empty.rst")
end

@testset "Read Restart File - h2o" begin
    # Test H2O Restart File
    atom_names, atom_masses, atom_coords, atom_types = read_rst("data/test_h2o.rst")
    @test atom_names == ["o", "h", "h"]
    @test atom_masses == [15.9994, 1.00794, 1.00794]
    @test atom_coords isa Matrix{Float64}
    @test atom_coords == [
        -0.01294293652 -0.00022738404 -0.00072569662
        0.43396012267 0.59193103593 0.66981595463
        0.68120631717 -0.49583217141 -0.52277201300
    ]
    @test atom_types == [1, 1, 1]
end

@testset "Read Restart File - nh3" begin
    # Test NH3 Restart File
    atom_names, atom_masses, atom_coords, atom_types = read_rst("data/test_nh3.rst")
    @test atom_names == ["N", "H", "H", "H"]
    @test atom_masses == [14.0067, 1.00794, 1.00794, 1.00794]
    @test atom_coords isa Matrix{Float64}
    @test atom_coords == [
        -0.000000 0.000000 0.084892
        0.000000 0.961804 -0.198082
        -0.832947 -0.480902 -0.198082
        0.832947 -0.480902 -0.198082
    ]
    @test atom_types == [2, 2, 2, 2]
end

@testset "Read Hessian File Exceptions" begin
    # Test Hessian File Exceptions
    @test_throws ErrorException read_hessian("data/not_a_file.dat")
    @test_throws ErrorException read_hessian("data/empty_hessian.dat")
    @test_throws ErrorException read_hessian("data/test_nh3.rst")
    @test_throws ErrorException read_hessian("data/test_hessian2.dat")

end

@testset "Read Hessian File" begin
    # Test Hessian file read correctly
    hessian_h2o = read_hessian("data/test_hessian_h2o.dat")
    @test hessian_h2o isa Matrix{Float64}
    @test size(hessian_h2o) == (9, 9)
    @test read_hessian("data/test_hessian.dat") == [2.0 3.0; 1.0 -2.0]
end

@testset "Read Moldescriptor Exceptions" begin
    atom_names = ["o", "h", "h"]
    atom_types = [1, 1, 1]
    # Test Moldescriptor File Exceptions
    @test_throws ErrorException read_moldescriptor("data/not_a_file.dat", atom_names, atom_types)
    @test_throws ErrorException read_moldescriptor("data/empty_moldescriptor.dat", atom_names, atom_types)
    @test_throws ErrorException read_moldescriptor("data/test_h2o.rst", atom_names, atom_types)
    @test_throws ErrorException read_moldescriptor("data/test_hessian.dat", atom_names, atom_types)
end

@testset "Read Moldescriptor File - h2o" begin
    # Test Moldescriptor file read correctly
    atom_names = ["o", "h", "h"]
    atom_types = [1, 1, 1]
    atom_charges = read_moldescriptor("data/test_moldescriptor.dat", atom_names, atom_types)
    @test atom_charges isa Vector{Float64}
    @test size(atom_charges) == (3,)
    @test atom_charges == [-0.65966, 0.32983, 0.32983]
end

@testset "Read Moldescriptor File - nh3" begin
    # Test Moldescriptor file read correctly
    atom_names = ["N", "H", "H", "H"]
    atom_types = [2, 2, 2, 2]
    atom_charges = read_moldescriptor("data/test_moldescriptor.dat", atom_names, atom_types)
    @test atom_charges isa Vector{Float64}
    @test size(atom_charges) == (4,)
    @test atom_charges == [-0.8022, 0.2674, 0.2674, 0.2674]
end