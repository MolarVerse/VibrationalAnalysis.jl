@testset "Read Files and Calculate VibAnal" begin
	atom_names, atom_masses, atom_coords, atom_types = read_rst("data/test_h2o.rst")
	hessian = read_hessian("data/test_hessian_h2o.dat")
	atom_charges = read_moldescriptor("data/test_moldescriptor.dat", atom_names, atom_types)
	wavenumbers, intensities, force_constants, reduced_masses = calculate(atom_masses, atom_coords, atom_charges, hessian)
	@test wavenumbers ≈ [0.134821, 0.173096, 0.190156, 36.56957, 41.79827, 49.84504, 1492.510456, 3669.03923, 3784.3906] atol = 1e-3
	@test intensities ≈ [17.66373, 17.27004, 17.54798, 143.56445, 94.1316, 168.41857, 196.35896, 114.02282, 212.83423] atol = 1e-3
	@test force_constants ≈ [0.0, 0.0, 0.0, 0.0, 0.0010716286, 0.0015643, 1.4173086, 8.31255, 9.171054] atol = 1e-3
	@test reduced_masses ≈ [6.00527, 6.005463, 6.0052792, 1.05927, 1.04103, 1.0686666, 1.079864, 1.0480, 1.0868] atol = 1e-4
end

@testset "Read Files and Calculate VibAnal - No Intensities" begin
	atom_names, atom_masses, atom_coords, atom_types = read_rst("data/test_h2o.rst")
	hessian = read_hessian("data/test_hessian_h2o.dat")
	wavenumbers, force_constants, reduced_masses = calculate(atom_masses, atom_coords, hessian)
	@test wavenumbers ≈ [0.134821, 0.173096, 0.190156, 36.56957, 41.79827, 49.84504, 1492.510456, 3669.03923, 3784.3906] atol = 1e-3
	@test force_constants ≈ [0.0, 0.0, 0.0, 0.0, 0.0010716286, 0.0015643, 1.4173086, 8.31255, 9.171054] atol = 1e-3
	@test reduced_masses ≈ [6.00527, 6.005463, 6.0052792, 1.05927, 1.04103, 1.0686666, 1.079864, 1.0480, 1.0868] atol = 1e-4
end

@testset "Normal Modes - H2O" begin
	atom_names, atom_masses, atom_coords, atom_types = read_rst("data/h2o_compare/final_h2o.xyz.rst")
	hessian = read_hessian("data/h2o_compare/hessian-5p-0.001.dat")
	_, _, _, eigenvectors_internal_normalized = calculate(atom_masses, atom_coords, hessian)
	test_modes = hcat(map(row -> parse.(Float64, row), split.(readlines("data/h2o_compare/normal_modes_h2o.dat")))...)
	@test eigenvectors_internal_normalized ≈ test_modes atol = 1e-5
end

@testset "Normal Modes - MeOH" begin
	atom_names, atom_masses, atom_coords, atom_types = read_rst("data/meoh_compare/final_meoh.xyz.rst")
	hessian = read_hessian("data/meoh_compare/hessian.dat")
	_, _, _, eigenvectors_internal_normalized = calculate(atom_masses, atom_coords, hessian)
	test_modes = hcat(map(row -> parse.(Float64, row), split.(readlines("data/meoh_compare/meoh_normal_modes.dat")))...)
	@test eigenvectors_internal_normalized ≈ test_modes atol = 1e-5
end
