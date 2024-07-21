using VibrationalAnalysis: wavenumber_hartree, wavenumber_eV, wavenumber_kcal, reduced_mass, force_constant
using Test

@testset "Wave Number hartree" begin
    wavenumber, omega = wavenumber_hartree([1.0])
    @test wavenumber ≈ [5140.4873] atol=1e-4
    @test omega ≈ [9.6829e14] rtol=1e-4
end

@testset "Wave Number eV" begin
    wavenumber, omega = wavenumber_eV([1.0])
    @test wavenumber ≈ [521.470831] atol=1e-4
    @test omega ≈ [9.82278e13] rtol=1e-4
end

@testset "Wave Number kcal" begin
    wavenumber, omega = wavenumber_kcal([1.0])
    @test wavenumber ≈ [108.5914] atol=1e-4
    @test omega ≈ [2.0455e13] rtol=1e-4 
end

@testset "Reduced Mass" begin
    @test reduced_mass([1.0, 2.0, 3.0]) == [1.0, 4.0, 9.0]
end

@testset "Force Constant" begin
    @test force_constant([1.0], [1.0]) == [1.6605778811026238e-29]

end

