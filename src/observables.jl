
export wavenumber_kcal, wavenumber_dftb, reduced_mass, force_constant, infrared_intensity

"""
    wavenumber_dftb(eigenvalues::Vector{Float64})

Convert eigenvalues from Hartree Å^-2 g^-1 to wavenumbers in cm^-1.

"""

function wavenumber_dftb(eigenvalues::Vector{Float64})
    
    # Conversion Hartree B^-2 g^-1 to s^-2: 
    #   2625500.2 J/Hartree * (0.188972598857892E+11)^2 B^2 / m^2  * 1000 g / kg
    omega = sqrt.(eigenvalues * 2625500.2 * (0.188972598857892E+11)^2)

    # Convert wavenumbers in s^-2 to wavenumbers in cm^-1 : 1/(2π * c) * ω
    wavenumbers = 1/(2π * 299792458) * omega

    return wavenumbers
end


"""
    wavenumber_kcal(eigenvalues::Vector{Float64})

Convert eigenvalues from kcal Å^-2 g^-1 to wavenumbers in cm^-1.

"""

function wavenumber_kcal(eigenvalues::Vector{Float64})
    
    # Convert eigenvalues in kcal Å^-2 g^-1 to wavenumbers in s^-2: 
    #   4184 * 10^23 (4184 J/kcal * 10^20 Å^2 / m^2 * 1000 g / kg)
    omega = sqrt.(eigenvalues * 4184.0 * 10e23)

    # Convert wavenumbers in s^-2 to wavenumbers in cm^-1 : 1/(2π * c) * ω
    wavenumbers = 1/(2π * 299792458) * omega

    return wavenumbers
end

"""
    reduced_mass(normalization::Vector{Float64})

Calculate the reduced mass from the normalization vector.
"""

function reduced_mass(normalization::Vector{Float64})
    return normalization.^2
end

"""
    force_constant(wavenumbers::Vector{Float64}, reduced_mass::Vector{Float64})

Calculate the force constant from the wavenumbers and the reduced mass.
"""

function force_constant(wavenumbers::Vector{Float64}, reduced_mass::Vector{Float64})
    return wavenumbers.^2 .* reduced_mass' / 6.022 / 1E28
end

"""
    infrared_intensity(normalization::Vector{Float64}, coordinates::Matrix{Float64}, charges::Vector{Float64})

Calculate the infrared intensity from the normalization vector, the coordinates and the charges.
"""

function infrared_intensity(normalization::Vector{Float64}, coordinates::Matrix{Float64}, charges::Vector{Float64})
    intensities = []
    for i in 1:size(normalization)[1]
        eigenvector = reshape(normalization[:, i], 3, :)'
        mu1 = sum((coordinates .+ 0.05eigenvector) .* charges, dims=1)
        mu2 = sum((coordinates .- 0.05eigenvector) .* charges, dims=1)
        delta_mu = mu1 - mu2 # norm
        intensity = sum((delta_mu ./ 2norm(0.1eigenvector)).^2)
        push!(intensities, intensity)
    end
    return intensities
end


