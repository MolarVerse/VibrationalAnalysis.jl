using Test
using VibrationalAnalysis

# Test reading RST files src/read_rst.jl

@testset "Read Restart File Exceptions" begin
    @test_throws ErrorException read_rst("../data/not_a_file.rst") 
    @test_throws ErrorException read_rst("../data/empty.rst")
end

@testset "Read Restart File" begin
    
    # Test H2O Restart File
    atom_names, atom_masses, atom_coords = read_rst("../data/test_h2o.rst")
    @test atom_names == ["o", "h", "h"]
    @test atom_masses == [15.9994, 1.00794, 1.00794]
    @test atom_coords isa Matrix{Float64}
    @test atom_coords == [ -0.01294293652 -0.00022738404 -0.00072569662;
                            0.43396012267  0.59193103593  0.66981595463;
                            0.68120631717 -0.49583217141 -0.52277201300 ] 

    # Test NH3 Restart File
    atom_names, atom_masses, atom_coords = read_rst("../data/test_nh3.rst")
    @test atom_names == ["N", "H", "H", "H"]
    @test atom_masses == [14.0067, 1.00794, 1.00794, 1.00794]
    @test atom_coords == [ -0.000000    0.000000    0.084892;
                            0.000000    0.961804   -0.198082;
                           -0.832947   -0.480902   -0.198082;
                            0.832947   -0.480902   -0.198082 ]
end


