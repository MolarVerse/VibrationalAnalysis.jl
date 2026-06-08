using VibrationalAnalysis: translational_modes, rotational_modes, transformation_matrix, internal_subspace, internal_coordinates
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

@testset "Linear molecule rotations" begin
	atom_coords = [0.0 0.0 -0.37; 0.0 0.0 0.37]
	atom_masses = [1.0, 1.0]
	rotation = rotational_modes(atom_coords, atom_masses)
	transformation = transformation_matrix(atom_coords, atom_masses)
	subspace = internal_subspace(atom_coords, atom_masses)
	@test size(rotation) == (6, 2)
	@test size(transformation) == (6, 5)
	@test size(subspace) == (6, 1)
	@test all(isfinite, rotation)
	@test all(isfinite, transformation)
	@test all(isfinite, subspace)
end

@testset "CO2 linear and near-linear rotations" begin
	atom_masses = [15.9994, 12.0107, 15.9994]

	linear_coords = [-1.16 0.0 0.0; 0.0 0.0 0.0; 1.16 0.0 0.0]
	for atom_coords in (
		linear_coords,
		[-1.16 0.0 0.0; 0.0 0.0 1e-6; 1.16 0.0 0.0],
		[-1.16 0.0 1e-6; 0.0 0.0 0.0; 1.16 0.0 0.0],
	)
		rotation = rotational_modes(atom_coords, atom_masses)
		transformation = transformation_matrix(atom_coords, atom_masses)
		subspace = internal_subspace(atom_coords, atom_masses)
		@test size(rotation) == (9, 2)
		@test size(transformation) == (9, 5)
		@test size(subspace) == (9, 4)
		@test all(isfinite, rotation)
		@test all(isfinite, transformation)
		@test all(isfinite, subspace)
	end

	bent_coords = [-1.16 0.0 0.0; 0.0 0.0 1e-5; 1.16 0.0 0.0]
	@test size(rotational_modes(bent_coords, atom_masses)) == (9, 3)
	@test size(transformation_matrix(bent_coords, atom_masses)) == (9, 6)
	@test size(internal_subspace(bent_coords, atom_masses)) == (9, 3)
end

@testset "Rotational Modes" begin
	atom_coords = [0.0 0.0 0.0; 0.0 0.0 1.0; 0.0 1.0 0.0]
	atom_masses = [1.0, 1.0, 1.0]
	rotation = rotational_modes(atom_coords, atom_masses)
	@test size(rotation) == (9, 3)
	@test rotation == [
		0.0 0.0 -0.2886751345948128;
		0.0 0.0 -0.2886751345948128;
		0.0 0.0 0.5773502691896257;
		-0.408248290463863 0.408248290463863 0.0;
		0.8164965809277261 0.408248290463863 0.0;
		-0.408248290463863 -0.8164965809277261 0.0;
		0.0 0.0 0.2886751345948128;
		0.0 0.0 -0.5773502691896257;
		0.0 0.0 0.2886751345948128
	]
end

@testset "Transformational Matrix" begin
	atom_coords = [0.0 0.0 0.0; 0.0 0.0 1.0; 0.0 1.0 0.0]
	atom_masses = [1.0, 1.0, 1.0]
	translation = translational_modes(atom_masses)
	rotation = rotational_modes(atom_coords, atom_masses)
	transformation = transformation_matrix(atom_coords, atom_masses)
	@test size(transformation) == (9, 6)
	@test transformation[:, 1:3] == translation
	@test transformation[:, 4:6] == rotation
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
