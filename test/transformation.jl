using VibrationalAnalysis: translational_modes, rotational_modes, transformation_matrix, internal_coordinates
using Test
using LinearAlgebra

@testset "Translation Modes" begin
	atom_masses = [1.0, 1.0, 1.0]
	translation = translational_modes(atom_masses)
	@test size(translation) == (9, 3)
	@test translation == [
		0.5773502691896258 0.0 0.0
		0.0 0.5773502691896258 0.0
		0.0 0.0 0.5773502691896258
		0.5773502691896258 0.0 0.0
		0.0 0.5773502691896258 0.0
		0.0 0.0 0.5773502691896258
		0.5773502691896258 0.0 0.0
		0.0 0.5773502691896258 0.0
		0.0 0.0 0.5773502691896258
	]
end

@testset "Rotational Modes" begin
	atom_coords = [0.0 0.0 0.0; 0.0 0.0 1.0; 0.0 1.0 0.0]
	atom_masses = [1.0, 1.0, 1.0]
	rotation = rotational_modes(atom_coords, atom_masses)
	@test size(rotation) == (9, 3)
	@test rotation == [
		-0.5773502691896257 0.0 0.0
		0.28867513459481303 -0.5 -0.43301270189221935
		0.28867513459481303 0.5 0.43301270189221935
		-0.5773502691896257 0.0 0.0
		0.28867513459481303 -0.5 0.43301270189221935
		0.28867513459481303 0.5 -0.43301270189221935
		0.0 0.0 0.40824829046386296
		0.0 0.0 -0.20412414523193156
		0.0 0.0 -0.20412414523193156
	]
end

@testset "Transformational Matrix" begin
	atom_coords = [0.0 0.0 0.0; 0.0 0.0 1.0; 0.0 1.0 0.0]
	atom_masses = [1.0, 1.0, 1.0]
	transformation = transformation_matrix(atom_coords, atom_masses)
	@test size(transformation) == (9, 9)
	@test transformation[:, 7:9] == zeros(9, 3)
end

@testset "Internal Coordinates" begin
	atom_coords = [0.0 0.0 0.0; 0.0 0.0 1.0; 0.0 1.0 0.0]
	atom_masses = [1.0, 1.0, 1.0]
	hessian_mw = diagm(0 => ones(9))
	eigenvalues, eigenvectors, normalization = internal_coordinates(atom_coords, atom_masses, hessian_mw)
	@test eigenvalues ≈ ones(9)
	@test size(eigenvectors) == (9, 9)
    @test normalization[:] ≈ ones(9)
end
