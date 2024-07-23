@testset "Write modes" begin
    temp = mktempdir()

    cd(temp) do
        eigenvectors_internal_normalized = hcat([1.0; 0.0; 0.0; 0.0; 1.0; 0.0])
        atom_coords = [1.0 0.0 0.0; 0.0 1.0 0.0]
        atom_names = ["H", "O"]
        write_modes(eigenvectors_internal_normalized, atom_coords, atom_names, filename="test_modes", amplitude=0.1, step=0.1)
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

@testset "Write wavenumber and intensities" begin
    wavenumbers = [1.0, 2.0, 3.0]
    intensities = [4.0, 5.0, 6.0]

    test_string = [
        "# Wavenumbers (cm-1)\t\tIntensities (km mol-1)",
        "1.0\t\t4.0",
        "2.0\t\t5.0",
        "3.0\t\t6.0"
    ]
    
    temp = mktempdir()
    cd(temp) do
        write_calculate_output(wavenumbers, intensities, filename="test_wavenumbers.dat")
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
        "# Wavenumbers (cm-1)\t\tIntensities (km mol-1)\t\tForce constants (mdyn Å-1)\t\tReduced masses (amu)",
        "1.0\t\t4.0\t\t7.0\t\t10.0",
        "2.0\t\t5.0\t\t8.0\t\t11.0",
        "3.0\t\t6.0\t\t9.0\t\t12.0"
    ]
    
    temp = mktempdir()
    cd(temp) do
        write_calculate_output(wavenumbers, intensities, force_constants, reduced_masses, filename="test_wavenumbers.dat")
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
        "# Wavenumbers (cm-1)\t\tForce constants (mdyn Å-1)\t\tReduced masses (amu)",
        "1.0\t\t7.0\t\t10.0",
        "2.0\t\t8.0\t\t11.0",
        "3.0\t\t9.0\t\t12.0"
    ]
    
    temp = mktempdir()
    cd(temp) do
        write_calculate_output(wavenumbers, force_constants, reduced_masses, filename="test_wavenumbers.dat")
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
