using VibrationalAnalysis
using Test

@testset "Translation Modes" begin
    atom_masses = [1.0]
    translation = translational_modes(atom_masses)
    @test size(translation) == (3, 3)
    @test translation == [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
end

@testset "Rotational Modes" begin
    atom_coords = [0.0 0.0 0.0; 0.0 0.0 1.0; 0.0 1.0 0.0]
    atom_masses = [1.0, 1.0, 1.0]
    rotation = rotational_modes(atom_coords, atom_masses)
    @test size(rotation) == (9, 3)
    @test rotation == 
                    [  -0.5773502691896257 0.0 0.0; 
                        0.28867513459481303 -0.5 -0.43301270189221935; 
                        0.28867513459481303 0.5 0.43301270189221935; 
                       -0.5773502691896257 0.0 0.0; 
                        0.28867513459481303 -0.5 0.43301270189221935; 
                        0.28867513459481303 0.5 -0.43301270189221935; 
                        0.0 0.0 0.40824829046386296; 
                        0.0 0.0 -0.20412414523193156; 
                        0.0 0.0 -0.20412414523193156
                    ]
end