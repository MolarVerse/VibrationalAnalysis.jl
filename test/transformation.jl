using VibrationalAnalysis
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
	@test transformation[:, 1] == [
		0.5773502691896258, 0.0, 0.0,
		0.5773502691896258, 0.0, 0.0,
		0.5773502691896258, 0.0, 0.0,
	]
	@test transformation[:, 2] == [
		0.0, 0.5773502691896258, 0.0,
		0.0, 0.5773502691896258, 0.0,
		0.0, 0.5773502691896258, 0.0,
	]
	@test transformation[:, 3] == [
		0.0, 0.0, 0.5773502691896258,
		0.0, 0.0, 0.5773502691896258,
		0.0, 0.0, 0.5773502691896258,
	]
	@test transformation[:, 4] == [
		-0.5773502691896257, 0.28867513459481303, 0.28867513459481303,
		-0.5773502691896257, 0.28867513459481303, 0.28867513459481303,
		0.0, 0.0, 0.0,
	]
	@test transformation[:, 5] == [
		0.0, -0.5, 0.5,
		0.0, -0.5, 0.5,
		0.0, 0.0, 0.0,
	]
	@test transformation[:, 6] == [
		0.0, -0.43301270189221935, 0.43301270189221935,
		0.0, 0.43301270189221935, -0.43301270189221935,
		0.40824829046386296, -0.20412414523193156, -0.20412414523193156,
	]
	@test transformation[:, 7:9] == zeros(9, 3)
end

@testset "Internal Coordinates" begin
	atom_coords = [0.0 0.0 0.0; 0.0 0.0 1.0; 0.0 1.0 0.0]
	atom_masses = [1.0, 1.0, 1.0]
	hessian_mw = diagm(0 => ones(9))
	eigenvalues, eigenvectors, normalization = internal_coordinates(atom_coords, atom_masses, hessian_mw)
	@test eigenvalues ≈ ones(9)
	@test eigenvectors[:, 1] == [
		-0.2787094021560701, -0.02758618133308042, 0.36021534146852785,
		0.12932248782060032, -0.20932598204740102, -0.595702806839151,
		0.1396447424776136, -0.5972568965401952, -0.013602141789088964,
	]
    @test eigenvectors[:, 9] == [
        0.4131993349694564, 0.33628890263202604, 0.023474557652313096, 
        0.20177621195797016, 0.08176280114106561, 0.19152455779067581, 
        -0.5796663737511388, -0.19199988522156133, 0.5085936652260151
    ]
    @test normalization[:] ≈ ones(9)
end
