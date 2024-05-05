@testset "Write wavenumber and frequency" begin
    temp = mktempdir()

    cd(temp) do
        wavenumbers = [1.0, 2.0, 3.0]
        intensities = [4.0, 5.0, 6.0]
        write_wavenumber_intensity(wavenumbers, intensities, filename="test_wavenumbers.dat")
        wavenumbers_test = readlines("test_wavenumbers.dat")
        @test wavenumbers_test[1] == "# Wavenumbers (cm-1)    Intensities (km mol-1)"
        @test wavenumbers_test[2] == "1.0    4.0"
        @test wavenumbers_test[3] == "2.0    5.0"
        @test wavenumbers_test[4] == "3.0    6.0"
    end
end

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

