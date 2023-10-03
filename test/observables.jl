using VibrationalAnalysis
using Test

@testset "Wave Number DFTB" begin
    @test wavenumber_dftb([1.0]) == [16255.648219071298]
end

@testset "Wave Number kcal" begin
    @test wavenumber_kcal([1.0]) == [34339.602736400884]
end

@testset "Reduced Mass" begin
    @test reduced_mass([1.0, 2.0, 3.0]) == [1.0, 4.0, 9.0]
end

@testset "Force Constant" begin
    @test force_constant([1.0], [1.0]) == [1.6605778811026238e-29]
end

