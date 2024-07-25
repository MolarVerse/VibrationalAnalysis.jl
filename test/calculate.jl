@testset "Read Files and Calculate VibAnal" begin
	atom_names, atom_masses, atom_coords, atom_types = read_rst("data/test_h2o.rst")
	hessian = read_hessian("data/test_hessian_h2o.dat")
	atom_charges = read_moldescriptor("data/test_moldescriptor.dat", atom_names, atom_types)
	wavenumbers, intensities, force_constants, reduced_masses = calculate(atom_masses, atom_coords, atom_charges, hessian)
	@test wavenumbers[7:end] ≈ [1492.5104569439752, 3669.03923547886, 3784.390612155495] atol = 1e-5
	@test intensities[7:end] ≈ [150.52494226058752, 86.42940396506292, 164.0700800418807] atol = 1e-5
	@test force_constants[7:end] ≈ [1.4173086016921352, 8.312559679465139, 9.171054294242243] atol = 1e-5
	@test reduced_masses[7:end] ≈ [1.0798645051916647, 1.0480202469197442, 1.086843369500246] atol = 1e-5
end

@testset "Read Files and Calculate VibAnal - No Intensities" begin
	atom_names, atom_masses, atom_coords, atom_types = read_rst("data/test_h2o.rst")
	hessian = read_hessian("data/test_hessian_h2o.dat")
	wavenumbers, force_constants, reduced_masses = calculate(atom_masses, atom_coords, hessian)
	@test wavenumbers[7:end] ≈ [1492.5104569439752, 3669.03923547886, 3784.390612155495] atol = 1e-5
	@test force_constants[7:end] ≈ [1.4173086016921352, 8.312559679465139, 9.171054294242243] atol = 1e-5
	@test reduced_masses[7:end] ≈ [1.0798645051916647, 1.0480202469197442, 1.086843369500246] atol = 1e-5
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
	@test abs.(eigenvectors_internal_normalized[:, 7:end]) ≈ abs.(test_modes[:, 7:end]) atol = 1e-5
end