import VibrationalAnalysis.check_unit

@testset "Unit Check" begin
    @test check_unit("kcal") == wavenumber_kcal
    @test check_unit("hartree") == wavenumber_hartree
    @test check_unit("ev") == wavenumber_eV
    @test check_unit("KcAl") == wavenumber_kcal
    @test check_unit("Hartree") == wavenumber_hartree
    @test check_unit("eV") == wavenumber_eV

    @test_throws ArgumentError check_unit("invalid")
end