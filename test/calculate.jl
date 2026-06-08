using LinearAlgebra
using VibrationalAnalysis: hessian_sign_factor, internal_subspace

@testset "Read Files and Calculate VibAnal" begin
	atom_names, atom_masses, atom_coords, atom_types = read_rst("data/test_h2o.rst")
	hessian = read_hessian("data/test_hessian_h2o.dat")
	atom_charges = read_moldescriptor("data/test_moldescriptor.dat", atom_names, atom_types)
	wavenumbers, intensities, force_constants, reduced_masses = calculate(atom_masses, atom_coords, atom_charges, hessian)
	@test wavenumbers[7:end] ≈ [1492.5348853402552, 3669.024005934286, 3784.327469724926] atol = 1e-5 rtol = 1e-8
	@test intensities[7:end] ≈ [150.5208654923484, 86.4334793788946, 164.08630456843233] atol = 1e-5 rtol = 1e-8
	@test force_constants[7:end] ≈ [1.417347274653204, 8.312499126458267, 9.170834256960058] atol = 1e-5 rtol = 1e-8
	@test reduced_masses[7:end] ≈ [1.079858621514245, 1.0480213128861526, 1.086853561248339] atol = 1e-5 rtol = 1e-8
end

@testset "Read Files and Calculate VibAnal - No Intensities" begin
	atom_names, atom_masses, atom_coords, atom_types = read_rst("data/test_h2o.rst")
	hessian = read_hessian("data/test_hessian_h2o.dat")
	wavenumbers, force_constants, reduced_masses = calculate(atom_masses, atom_coords, hessian)
	@test wavenumbers[7:end] ≈ [1492.5348853402552, 3669.024005934286, 3784.327469724926] atol = 1e-5 rtol = 1e-8
	@test force_constants[7:end] ≈ [1.417347274653204, 8.312499126458267, 9.170834256960058] atol = 1e-5 rtol = 1e-8
	@test reduced_masses[7:end] ≈ [1.079858621514245, 1.0480213128861526, 1.086853561248339] atol = 1e-5 rtol = 1e-8
end

@testset "Linear molecule calculation" begin
	atom_masses = [1.0, 1.0]
	atom_coords = [0.0 0.0 -0.37; 0.0 0.0 0.37]
	hessian = Matrix{Float64}(I, 6, 6)
	wavenumbers, force_constants, reduced_masses, modes = calculate(atom_masses, atom_coords, hessian, hessian_sign = :positive)
	@test size(modes) == (6, 6)
	@test length(wavenumbers) == 6
	@test length(force_constants) == 6
	@test length(reduced_masses) == 6
	@test all(isfinite, wavenumbers)
	@test all(isfinite, force_constants)
	@test all(isfinite, reduced_masses)
	@test all(isfinite, modes)
end

@testset "CO2 linear and near-linear calculation" begin
	atom_masses = [15.9994, 12.0107, 15.9994]
	hessian = Matrix{Float64}(I, 9, 9)

	for atom_coords in (
		[-1.16 0.0 0.0; 0.0 0.0 0.0; 1.16 0.0 0.0],
		[-1.16 0.0 0.0; 0.0 0.0 1e-6; 1.16 0.0 0.0],
		[-1.16 0.0 1e-6; 0.0 0.0 0.0; 1.16 0.0 0.0],
		[-1.16 0.0 0.0; 0.0 0.0 1e-5; 1.16 0.0 0.0],
	)
		wavenumbers, force_constants, reduced_masses, modes = calculate(atom_masses, atom_coords, hessian, hessian_sign = :positive)
		@test size(modes) == (9, 9)
		@test length(wavenumbers) == 9
		@test length(force_constants) == 9
		@test length(reduced_masses) == 9
		@test all(isfinite, wavenumbers)
		@test all(isfinite, force_constants)
		@test all(isfinite, reduced_masses)
		@test all(isfinite, modes)
	end
end

@testset "Hessian sign convention" begin
	atom_masses = [1.0, 1.0]
	atom_coords = [0.0 0.0 -0.37; 0.0 0.0 0.37]
	hessian = -Matrix{Float64}(I, 6, 6)
	auto_wavenumbers, _, _, _ = calculate(atom_masses, atom_coords, hessian)
	positive_wavenumbers, _, _, _ = calculate(atom_masses, atom_coords, hessian, hessian_sign = :positive)
	negative_wavenumbers, _, _, _ = calculate(atom_masses, atom_coords, -hessian, hessian_sign = :negative)
	numeric_wavenumbers, _, _, _ = calculate(atom_masses, atom_coords, hessian, hessian_sign = -1)
	@test all(isfinite, auto_wavenumbers)
	@test all(isfinite, positive_wavenumbers)
	@test all(isfinite, negative_wavenumbers)
	@test all(isfinite, numeric_wavenumbers)
	@test all(>=(0), auto_wavenumbers)
	@test any(<(0), positive_wavenumbers)
	@test any(<(0), negative_wavenumbers)
	@test all(>=(0), numeric_wavenumbers)
	@test_throws ArgumentError calculate(atom_masses, atom_coords, hessian, hessian_sign = :invalid)

	single_atom_wavenumbers, _, _, single_atom_modes = calculate([1.0], [0.0 0.0 0.0], Matrix{Float64}(I, 3, 3))
	@test size(single_atom_modes) == (3, 3)
	@test all(isfinite, single_atom_wavenumbers)

	triatomic_coords = [0.0 0.0 0.0; 0.0 0.0 1.0; 0.0 1.0 0.0]
	triatomic_masses = [1.0, 1.0, 1.0]
	subspace = internal_subspace(triatomic_coords, triatomic_masses)
	tie_hessian = subspace * Diagonal([-2.0, 1.0, 0.0]) * subspace'
	@test hessian_sign_factor(triatomic_coords, triatomic_masses, tie_hessian, :auto) == -1.0
end

@testset "Normal Modes - H2O" begin
	atom_names, atom_masses, atom_coords, atom_types = read_rst("data/h2o_compare/final_h2o.xyz.rst")
	hessian = read_hessian("data/h2o_compare/hessian-5p-0.001.dat")
	_, _, _, eigenvectors_internal_normalized = calculate(atom_masses, atom_coords, hessian)
	test_modes = hcat(map(row -> parse.(Float64, row), split.(readlines("data/h2o_compare/normal_modes_h2o.dat")))...)'
	@test abs.(eigenvectors_internal_normalized[:, 7:end]) ≈ abs.(test_modes[:, 7:end]) atol = 1e-5
end

@testset "Normal Modes - MeOH" begin
	atom_names, atom_masses, atom_coords, atom_types = read_rst("data/meoh_compare/final_meoh.xyz.rst")
	hessian = read_hessian("data/meoh_compare/hessian.dat")
	_, _, _, eigenvectors_internal_normalized = calculate(atom_masses, atom_coords, hessian)
	test_modes = hcat(map(row -> parse.(Float64, row), split.(readlines("data/meoh_compare/meoh_normal_modes.dat")))...)'
	mode_diffs = abs.(abs.(eigenvectors_internal_normalized[:, 7:end]) .- abs.(test_modes[:, 7:end]))
	@test maximum(mode_diffs) <= 1.1e-4
end
