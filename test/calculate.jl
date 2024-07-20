@testset "Read Files and Calculate VibAnal" begin
	rst_file = "data/test_h2o.rst"
	hessian_file = "data/test_hessian_h2o.dat"
	moldescriptor_file = "data/test_moldescriptor.dat"
	wavenumbers, intensities, force_constants, reduced_masses = read_calculate(rst_file, hessian_file, moldescriptor_file)
	@test wavenumbers ≈ [0.134821, 0.173096, 0.190156, 36.56957, 41.79827, 49.84504, 1492.510456, 3669.03923, 3784.3906] atol = 1e-3
	@test intensities ≈ [17.66373, 17.27004, 17.54798, 143.56445, 94.1316, 168.41857, 196.35896, 114.02282, 212.83423] atol = 1e-3
	@test force_constants ≈ [0.0, 0.0, 0.0, 0.0, 0.0010716286, 0.0015643, 1.4173086, 8.31255, 9.171054] atol = 1e-3
	@test reduced_masses ≈ [6.00527, 6.005463, 6.0052792, 1.05927, 1.04103, 1.0686666, 1.079864, 1.0480, 1.0868] atol = 1e-4
end

@testset "Read Files and Calculate VibAnal - No Intensities" begin
	rst_file = "data/test_h2o.rst"
	hessian_file = "data/test_hessian_h2o.dat"
	wavenumbers, intensities, force_constants, reduced_masses = read_calculate(rst_file, hessian_file)
	@test wavenumbers ≈ [0.134821, 0.173096, 0.190156, 36.56957, 41.79827, 49.84504, 1492.510456, 3669.03923, 3784.3906] atol = 1e-3
	@test force_constants ≈ [0.0, 0.0, 0.0, 0.0, 0.0010716286, 0.0015643, 1.4173086, 8.31255, 9.171054] atol = 1e-3
	@test reduced_masses ≈ [6.00527, 6.005463, 6.0052792, 1.05927, 1.04103, 1.0686666, 1.079864, 1.0480, 1.0868] atol = 1e-4
end
