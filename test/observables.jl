using VibrationalAnalysis: wavenumber_dftb, wavenumber_kcal, reduced_mass, force_constant
using Test

@testset "Wave Number DFTB" begin
    wavenumber, omega = wavenumber_dftb([1.0])
    @test wavenumber == [162.55648219071298]
    @test omega == [3.061997724547502e13]
end

@testset "Wave Number kcal" begin
    wavenumber, omega = wavenumber_kcal([1.0])
    @test wavenumber == [108.59135859237747]
    @test omega == [2.0454828280872953e13]
end

@testset "Reduced Mass" begin
    @test reduced_mass([1.0, 2.0, 3.0]) == [1.0, 4.0, 9.0]
end

@testset "Force Constant" begin
    @test force_constant([1.0], [1.0]) == [1.6605778811026238e-29]

end

