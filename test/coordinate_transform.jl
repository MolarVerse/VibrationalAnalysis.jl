using VibrationalAnalysis
using Test

@testset "Center of Mass - h2o" begin
    # Test H2O Restart File
    _ , atom_masses, atom_coords = read_rst("../data/test_h2o.rst")
    @test atom_masses == [15.9994, 1.00794, 1.00794]
    @test atom_coords == [
        -0.01294293652 -0.00022738404 -0.00072569662
        0.43396012267 0.59193103593 0.66981595463
        0.68120631717 -0.49583217141 -0.52277201300
    ]
    # Test center of mass
    @test center_to_com(atom_coords, atom_masses) == [
        -0.0638409321556072 -0.005402095801055762 -0.00830819103721218
        0.3830621270343928 0.5867563241689442 0.6622334602127878
        0.6303083215343929 -0.5010068831710558 -0.5303545074172121
    ]
end

@testset "Center of Mass - nh3" begin
    # Test NH3 Restart File
    _ , atom_masses, atom_coords = read_rst("../data/test_nh3.rst")
    @test atom_masses == [14.0067, 1.00794, 1.00794, 1.00794]
    @test atom_coords == [
        -0.000000 0.000000 0.084892
        0.000000 0.961804 -0.198082
        -0.832947 -0.480902 -0.198082
        0.832947 -0.480902 -0.198082
    ]
    # Test center of mass
    @test center_to_com(atom_coords, atom_masses) == [
        -0.0 0.0 0.05024288399179827
        0.0 0.961804 -0.23273111600820173
        -0.832947 -0.480902 -0.23273111600820173
        0.832947 -0.480902 -0.23273111600820173
    ]
end

@testset "Inertia Tensor - h2o" begin
    # Test H2O Restart File
    _ , atom_masses, atom_coords = read_rst("../data/test_h2o.rst")
    @test atom_masses == [15.9994, 1.00794, 1.00794]
    @test atom_coords == [
        -0.01294293652 -0.00022738404 -0.00072569662
        0.43396012267 0.59193103593 0.66981595463
        0.68120631717 -0.49583217141 -0.52277201300
    ]
    # Test inertia tensor
    @test inertia_tensor(atom_coords, atom_masses) == [
        1.3271332725386336 0.08622962761730685 0.07276422299611818;
        0.08622962761730685 1.3402017990291517 -0.6601939995475998;
        0.07276422299611818 -0.6601939995475998 1.2140373169840863
    ]
end

@testset "Inertia Tensor - nh3" begin
    # Test NH3 Restart File
    _ , atom_masses, atom_coords = read_rst("../data/test_nh3.rst")
    @test atom_masses == [14.0067, 1.00794, 1.00794, 1.00794]
    @test atom_coords == [
        -0.000000 0.000000 0.084892
        0.000000 0.961804 -0.198082
        -0.832947 -0.480902 -0.198082
        0.832947 -0.480902 -0.198082
    ]
    # Test inertia tensor
    @test inertia_tensor(atom_coords, atom_masses) == [
        1.5977572235586526 -0.0 -0.0;
        -0.0 1.5977582395561254 -0.0;
        -0.0 -0.0 2.7972369136232618
    ]
end

