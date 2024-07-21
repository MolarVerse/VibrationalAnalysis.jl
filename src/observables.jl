export wavenumber_hartree, wavenumber_eV, wavenumber_kcal

"""
	wavenumber_hartree(eigenvalues::Vector{Float64}) -> wavenumbers::Vector{Float64}, omega::Vector{Float64}

Convert `eigenvalues` from hartree bor^-2 g^-1 to `wavenumbers` in cm^-1. Made for DFTB hessian files.
Output include `omega` in s^-2.

``ω = √v``

``ν̃ = 1/(2π * c) * ω``

# Arguments
- `eigenvalues::Vector{Float64}`: The eigenvalues.
"""
function wavenumber_hartree(eigenvalues::Vector{Float64})

	# Conversion hartree bor^-2 g^-1 to s^-2: 
	#   2625500.2 J/Hartree * (0.188972598857892E+11)^2 bor^2 / m^2  * 1000 g / kg
	omega = sqrt.(eigenvalues * 2625500.2 * (0.188972598857892E+11)^2 * 1000)

	# Convert wavenumbers in s^-2 to wavenumbers in cm^-1 : 1/(2π * c) * ω
	wavenumbers = 1 / (2π * 2.99792458e10) * omega

	return wavenumbers, omega
end

"""
	wavenumber_eV(eigenvalues::Vector{Float64}) -> wavenumbers::Vector{Float64}, omega::Vector{Float64}

Convert `eigenvalues` from eV Å^-2 g^-1 to `wavenumbers` in cm^-1. Made for ASE hessian files.
Output include `omega` in s^-2.

``ω = √v``

``ν̃ = 1/(2π * c) * ω``

# Arguments
- `eigenvalues::Vector{Float64}`: The eigenvalues.
"""
function wavenumber_eV(eigenvalues::Vector{Float64})

	# Conversion eV Å^-2 g^-1 to s^-2: 
	#   96485.307499 J/eV * 10^20 Å^2 / m^2 * 1000 g / kg)
	omega = sqrt.(eigenvalues * 96485.307499 * 1.0E23)

	# Convert wavenumbers in s^-2 to wavenumbers in cm^-1 : 1/(2π * c) * ω
	wavenumbers = 1 / (2π * 2.99792458e10) * omega

	return wavenumbers, omega
end

"""
	wavenumber_kcal(eigenvalues::Vector{Float64}) -> wavenumbers::Vector{Float64}, omega::Vector{Float64}

Convert `eigenvalues` from kcal Å^-2 g^-1 to `wavenumbers` in cm^-1. Made for QMCFC hessian files.
Output include `omega` in s^-2.

``ω = √v``

``ν̃ = 1/(2π * c) * ω``

# Arguments
- `eigenvalues::Vector{Float64}`: The eigenvalues.
"""
function wavenumber_kcal(eigenvalues::Vector{Float64})

	# Convert eigenvalues in kcal Å^-2 g^-1 to wavenumbers in s^-2: 
	#   4184 * 10^23 (4184 J/kcal * 10^20 Å^2 / m^2 * 1000 g / kg)
	omega = sqrt.(eigenvalues * 4184.0 * 1.0E23)

	# Convert wavenumbers in s^-2 to wavenumbers in cm^-1 : 1/(2π * c) * ω
	wavenumbers = 1 / (2π * 2.99792458e10) * omega

	return wavenumbers, omega
end

"""
	reduced_mass(normalization::Vector{Float64}) -> red_mass::Vector{Float64}

Calculate the reduced mass from the normalization vector.

# Arguments
- `normalization::Vector{Float64}`: The normalization vector.
"""
function reduced_mass(normalization::Vector{Float64})
	red_mass = normalization .^ 2
	return red_mass[:] # convert to vector
end

"""
	force_constant(wavenumbers::Vector{Float64}, reduced_mass::Vector{Float64}) -> force_const::Vector{Float64}

Calculate the force constant from the wavenumbers and the reduced mass.

# Arguments
- `wavenumbers::Vector{Float64}`: The wavenumbers.
- `reduced_mass::Vector{Float64}`: The reduced mass.
"""
function force_constant(omega::Vector{Float64}, reduced_mass::Vector{Float64})
	# Conversion g mol^-1 s^-2 to mdyn Å^-1: / 6.022 / 1E23 (mol) / 1E3 (kg/g) / 1E2 (mdyn/Å / N/m)
	force_const = omega .^ 2 .* reduced_mass / 6.022 / 1E28
	return force_const[:] # convert to vector
end

"""
	infrared_intensity(eigenvectors_internal_normalized::Matrix{Float64}, atom_charges::Vector{Float64}, reduced_mass::Vector{Float64}) -> intensities::Vector{Float64}

Calculate the infrared intensity in km mol^-1 from the normalization eigen matrix, the coordinates and the charges.

# Arguments
- `eigenvectors_internal_normalized::Matrix{Float64}`: The normalized eigenvectors.
- `atom_charges::Vector{Float64}`: The atom charges.
- `reduced_mass::Vector{Float64}`: The reduced mass.
"""
function infrared_intensity(eigenvectors_internal_normalized::Matrix{Float64}, atom_charges::Vector{Float64}, reduced_masses::Vector{Float64})

	intensities = []

	for i in 1:size(eigenvectors_internal_normalized)[1]

		# Eigenvector of the i-th mode
		eigenvector = reshape(eigenvectors_internal_normalized[:, i], 3, :)'

		# Conversion factor taken from: http://thiele.ruc.dk/~spanget/help/g09/k_constants.html
		intensity = sum((sum(eigenvector .* atom_charges, dims = 1) / 0.2081943 / norm(eigenvector)) .^ 2 / reduced_masses[i] * 42.2561)

		push!(intensities, intensity)
	end

	return intensities
end
