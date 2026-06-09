using Test

@testset "CLI structure input - xyz" begin
    temp = mktempdir()
    xyz_file = joinpath(@__DIR__, "data", "test_h2o.xyz")
    hessian_file = joinpath(@__DIR__, "data", "test_hessian_h2o.dat")
    moldescriptor_file = joinpath(@__DIR__, "data", "test_h2o_moldescriptor.dat")

    cd(temp) do
        VibrationalAnalysis.vibrationalanalysis(xyz_file, hessian_file, output = "xyz_wavenumbers.dat")
        wavenumber_lines = readlines("xyz_wavenumbers.dat")
        @test wavenumber_lines[1] == "# Wavenumbers (cm-1)  Force constants (mdyn Å-1)  Reduced masses (amu)"
        @test length(wavenumber_lines) == 10

        VibrationalAnalysis.vibrationalanalysis(xyz_file, hessian_file, moldescriptor = moldescriptor_file, output = "xyz_intensities.dat")
        intensity_lines = readlines("xyz_intensities.dat")
        @test intensity_lines[1] == "# Wavenumbers (cm-1)  Intensities (km mol-1)  Force constants (mdyn Å-1)  Reduced masses (amu)"
        @test length(intensity_lines) == 10
    end
end

@testset "CLI structure input - xyz moldescriptor restrictions" begin
    xyz_file = joinpath(@__DIR__, "data", "test_h2o.xyz")
    hessian_file = joinpath(@__DIR__, "data", "test_hessian_h2o.dat")
    moldescriptor_file = joinpath(@__DIR__, "data", "test_moldescriptor.dat")
    @test_throws ErrorException VibrationalAnalysis.vibrationalanalysis(xyz_file, hessian_file, moldescriptor = moldescriptor_file)
end
