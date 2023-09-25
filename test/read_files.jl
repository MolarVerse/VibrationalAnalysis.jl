using VibrationalAnalysis
using Test

# Test reading RST files src/read_files.jl

@testset "Read Restart File Exceptions" begin
    @test_throws ErrorException read_rst("../data/not_a_file.rst")
    @test_throws ErrorException read_rst("../data/empty.rst")
end

@testset "Read Restart File - h2o" begin
    # Test H2O Restart File
    atom_names, atom_masses, atom_coords,  = read_rst("../data/test_h2o.rst")
    @test atom_names == ["o", "h", "h"]
    @test atom_masses == [15.9994, 1.00794, 1.00794]
    @test atom_coords isa Matrix{Float64}
    @test atom_coords == [
        -0.01294293652 -0.00022738404 -0.00072569662
        0.43396012267 0.59193103593 0.66981595463
        0.68120631717 -0.49583217141 -0.52277201300
    ]
end

@testset "Read Restart File - nh3" begin
    # Test NH3 Restart File
    atom_names, atom_masses, atom_coords,  = read_rst("../data/test_nh3.rst")
    @test atom_names == ["N", "H", "H", "H"]
    @test atom_masses == [14.0067, 1.00794, 1.00794, 1.00794]
    @test atom_coords == [
        -0.000000 0.000000 0.084892
        0.000000 0.961804 -0.198082
        -0.832947 -0.480902 -0.198082
        0.832947 -0.480902 -0.198082
    ]
end

@testset "Read Hessian File Exceptions" begin

    # Test Hessian File Exceptions
    @test_throws ErrorException read_hessian("../data/not_a_file.dat")
    @test_throws ErrorException read_hessian("../data/empty_hessian.dat")
    @test_throws ErrorException read_hessian("../data/test_nh3.rst")
    @test_throws ErrorException read_hessian("../data/test_hessian2.dat")

    # Test Hessian file read correctly
    hessian_h2o = read_hessian("../data/test_hessian_h2o.dat")
    @test hessian_h2o isa Matrix{Float64}
    @test size(hessian_h2o) == (9, 9)
    @test read_hessian("../data/test_hessian.dat") == [2.0 3.0; 1.0 -2.0]

end