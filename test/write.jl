@testset "Write modes" begin
	temp = mktempdir()

	cd(temp) do
		eigenvectors_internal_normalized = hcat([1.0; 0.0; 0.0; 0.0; 1.0; 0.0])
		atom_coords = [1.0 0.0 0.0; 0.0 1.0 0.0]
		atom_names = ["H", "O"]
		write_modes(eigenvectors_internal_normalized, atom_coords, atom_names, filename = "test_modes", amplitude = 0.1, step = 0.1)
		mode1 = readlines("test_modes-1.xyz")
		@test mode1[1] == "2"
		@test mode1[3] == "H    0.9   0.0   0.0"
		@test mode1[4] == "O    0.0   0.9   0.0"
		@test mode1[5] == "2"
		@test mode1[7] == "H    1.0   0.0   0.0"
		@test mode1[8] == "O    0.0   1.0   0.0"
		@test mode1[9] == "2"
		@test mode1[11] == "H    1.1   0.0   0.0"
		@test mode1[12] == "O    0.0   1.1   0.0"
	end
end

@testset "Write modes file" begin
	temp = mktempdir()

	cd(temp) do
		eigenvectors_internal_normalized = hcat([ 
                1.0 0.0 1.0 0.0 1.0 0.0;
				0.0 1.0 0.0 1.0 0.0 1.0;
				0.0 0.0 0.0 0.0 0.0 0.0;
				0.0 0.0 0.0 0.0 0.0 1.0;
				0.0 0.0 0.0 1.0 0.0 0.0;
				0.0 0.0 0.0 0.0 1.0 0.0
			])
		write_modes(eigenvectors_internal_normalized, filename = "modes.dat")
		mode1 = readlines("modes.dat")
		@test mode1[1] == "1.0 0.0 1.0 0.0 1.0 0.0 "
		@test mode1[2] == "0.0 1.0 0.0 1.0 0.0 1.0 "
		@test mode1[3] == "0.0 0.0 0.0 0.0 0.0 0.0 "
		@test mode1[4] == "0.0 0.0 0.0 0.0 0.0 1.0 "
		@test mode1[5] == "0.0 0.0 0.0 1.0 0.0 0.0 "
		@test mode1[6] == "0.0 0.0 0.0 0.0 1.0 0.0 "
	end
end

@testset "Write wavenumber and intensities" begin
	wavenumbers = [1.0, 2.0, 3.0]
	intensities = [4.0, 5.0, 6.0]

	test_string = [
		"# Wavenumbers (cm-1)  Intensities (km mol-1)",
		"1.00000000e+00\t4.00000000e+00",
		"2.00000000e+00\t5.00000000e+00",
		"3.00000000e+00\t6.00000000e+00",
	]

	temp = mktempdir()
	cd(temp) do
		write_calculate_output(wavenumbers, intensities, filename = "test_wavenumbers.dat")
		wavenumbers_test = readlines("test_wavenumbers.dat")
		@test wavenumbers_test[1] == test_string[1]
		@test wavenumbers_test[2] == test_string[2]
		@test wavenumbers_test[3] == test_string[3]
		@test wavenumbers_test[4] == test_string[4]

		open("stdout", "w") do io
			redirect_stdout(io) do
				write_calculate_output(wavenumbers, intensities)
			end
		end
		wavenumbers_test = readlines("stdout")
		@test wavenumbers_test[1] == test_string[1]
		@test wavenumbers_test[2] == test_string[2]
		@test wavenumbers_test[3] == test_string[3]
		@test wavenumbers_test[4] == test_string[4]
	end
end

@testset "Write calculate output - w/ I f μ" begin
	wavenumbers = [1.0, 2.0, 3.0]
	intensities = [4.0, 5.0, 6.0]
	force_constants = [7.0, 8.0, 9.0]
	reduced_masses = [10.0, 11.0, 12.0]

	test_string = [
		"# Wavenumbers (cm-1)  Intensities (km mol-1)  Force constants (mdyn Å-1)  Reduced masses (amu)",
		"1.00000000e+00\t4.00000000e+00\t7.00000000e+00\t1.00000000e+01",
		"2.00000000e+00\t5.00000000e+00\t8.00000000e+00\t1.10000000e+01",
		"3.00000000e+00\t6.00000000e+00\t9.00000000e+00\t1.20000000e+01",
	]

	temp = mktempdir()
	cd(temp) do
		write_calculate_output(wavenumbers, intensities, force_constants, reduced_masses, filename = "test_wavenumbers.dat")
		wavenumbers_test = readlines("test_wavenumbers.dat")
		@test wavenumbers_test[1] == test_string[1]
		@test wavenumbers_test[2] == test_string[2]
		@test wavenumbers_test[3] == test_string[3]
		@test wavenumbers_test[4] == test_string[4]

		open("stdout", "w") do io
			redirect_stdout(io) do
				write_calculate_output(wavenumbers, intensities, force_constants, reduced_masses)
			end
		end
		wavenumbers_test = readlines("stdout")
		@test wavenumbers_test[1] == test_string[1]
		@test wavenumbers_test[2] == test_string[2]
		@test wavenumbers_test[3] == test_string[3]
		@test wavenumbers_test[4] == test_string[4]
	end
end

@testset "Write calculate output - w/ f μ" begin
	wavenumbers = [1.0, 2.0, 3.0]
	force_constants = [7.0, 8.0, 9.0]
	reduced_masses = [10.0, 11.0, 12.0]

	test_string = [
		"# Wavenumbers (cm-1)  Force constants (mdyn Å-1)  Reduced masses (amu)",
		"1.00000000e+00\t7.00000000e+00\t1.00000000e+01",
		"2.00000000e+00\t8.00000000e+00\t1.10000000e+01",
		"3.00000000e+00\t9.00000000e+00\t1.20000000e+01",
	]

	temp = mktempdir()
	cd(temp) do
		write_calculate_output(wavenumbers, force_constants, reduced_masses, filename = "test_wavenumbers.dat")
		wavenumbers_test = readlines("test_wavenumbers.dat")
		@test wavenumbers_test[1] == test_string[1]
		@test wavenumbers_test[2] == test_string[2]
		@test wavenumbers_test[3] == test_string[3]
		@test wavenumbers_test[4] == test_string[4]

		open("stdout", "w") do io
			redirect_stdout(io) do
				write_calculate_output(wavenumbers, force_constants, reduced_masses)
			end
		end
		wavenumbers_test = readlines("stdout")
		@test wavenumbers_test[1] == test_string[1]
		@test wavenumbers_test[2] == test_string[2]
		@test wavenumbers_test[3] == test_string[3]
		@test wavenumbers_test[4] == test_string[4]
	end
end
