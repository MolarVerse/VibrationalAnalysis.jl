export calculate

"""
	calculate(atom_masses, atom_coords, hessian)

Calculates the wavenumbers, force constants and reduced masses from the atom masses, atom coordinates, and the hessian.

# Arguments
- `atom_masses::Vector{Float64}`: Vector of atom_masses (n)
- `atom_coords::Matrix{Float64}`: Matrix of atom coordinates (3xn)
- `hessian::Matrix{Float64}`: Matrix of the hessian (3nx3n)

# Keyword Arguments
- `wavenumber::Function`: The wavenumber function to use. Either `wavenumber_kcal`, `wavenumber_eV` or `wavenumber_hartree`. Default is `wavenumber_kcal`.
- `hessian_sign::Union{Symbol,Real}`: Hessian sign convention. Use `:auto`, `:positive`, `:negative`, `1`, or `-1`. Default is `:auto`.

# Returns
- `wavenumbers::Vector{Float64}`: The wavenumbers.
- `force_constants::Vector{Float64}`: The force constants.
- `reduced_masses::Vector{Float64}`: The reduced masses.
- `eigenvectors_internal_normalized::Matrix{Float64}`: The eigenvectors in internal coordinates.

# Example
```julia-repl
julia> calculate(atom_masses, atom_coords, hessian)
julia> calculate(atom_masses, atom_coords, hessian, wavenumber=wavenumber_hartree)
julia> calculate(atom_masses, atom_coords, hessian, wavenumber=wavenumber_eV)
```
"""
function calculate(atom_masses::Vector{Float64}, atom_coords::Matrix{Float64}, hessian::Matrix{Float64}; wavenumber = wavenumber_kcal, hessian_sign = :auto)

	# Symmetrize the hessian and mass weighting
	sign_factor = hessian_sign_factor(atom_coords, atom_masses, hessian, hessian_sign)
	hessian_mw = mass_weighted_hessian(hessian, atom_masses, sign = sign_factor)

	# Transform in internal coordinates
	eigenvalues, eigenvectors_internal_normalized, normalization = internal_coordinates(atom_coords, atom_masses, hessian_mw)

	# Calculate observables
	wavenumbers, omega = wavenumber(eigenvalues)
	reduced_masses = reduced_mass(normalization)
	force_constants = force_constant(omega, reduced_masses)

	return wavenumbers, force_constants, reduced_masses, eigenvectors_internal_normalized
end

"""
	calculate(atom_masses, atom_coords, atom_charges, hessian)

Calculates the wavenumbers, intensities, force constants and reduced masses from the atom masses, atom coordinates, atom charges and the hessian.

# Arguments
- `atom_masses::Vector{Float64}`: Vector of atom_masses (n)
- `atom_coords::Matrix{Float64}`: Matrix of atom coordinates (3xn)
- `atom_charges::Vector{Float64}`: Vector of atom charges (n)
- `hessian::Matrix{Float64}`: Matrix of the hessian (3nx3n)

# Keyword Arguments
- `wavenumber::Function`: The wavenumber function to use. Either `wavenumber_kcal`, `wavenumber_eV` or `wavenumber_hartree`. Default is `wavenumber_kcal`.
- `hessian_sign::Union{Symbol,Real}`: Hessian sign convention. Use `:auto`, `:positive`, `:negative`, `1`, or `-1`. Default is `:auto`.

# Returns
- `wavenumbers::Vector{Float64}`: The wavenumbers.
- `intensities::Vector{Float64}`: The intensities.
- `force_constants::Vector{Float64}`: The force constants.
- `reduced_masses::Vector{Float64}`: The reduced masses.
- `eigenvectors_internal_normalized::Matrix{Float64}`: The eigenvectors in internal coordinates.

# Example
```julia-repl
julia> calculate(atom_masses, atom_coords, atom_charges, hessian)
julia> calculate(atom_masses, atom_coords, atom_charges, hessian, wavenumber=wavenumber_hartree)
julia> calculate(atom_masses, atom_coords, atom_charges, hessian, wavenumber=wavenumber_eV)
```
"""
function calculate(atom_masses::Vector{Float64}, atom_coords::Matrix{Float64}, atom_charges::Vector{Float64}, hessian::Matrix{Float64}; wavenumber = wavenumber_kcal, hessian_sign = :auto)

	# Symmetrize the hessian and mass weighting
	sign_factor = hessian_sign_factor(atom_coords, atom_masses, hessian, hessian_sign)
	hessian_mw = mass_weighted_hessian(hessian, atom_masses, sign = sign_factor)

	# Transform in internal coordinates
	eigenvalues, eigenvectors_internal_normalized, normalization = internal_coordinates(atom_coords, atom_masses, hessian_mw)

	# Calculate observables
	wavenumbers, omega = wavenumber(eigenvalues)
	reduced_masses = reduced_mass(normalization)
	intensities = infrared_intensity(eigenvectors_internal_normalized, atom_charges, reduced_masses)
	force_constants = force_constant(omega, reduced_masses)

	return wavenumbers, intensities, force_constants, reduced_masses, eigenvectors_internal_normalized
end

function hessian_sign_factor(atom_coords::Matrix{Float64}, atom_masses::Vector{Float64}, hessian::Matrix{Float64}, hessian_sign)
	if hessian_sign isa Real
		if hessian_sign == 1 || hessian_sign == -1
			return Float64(hessian_sign)
		end
	elseif hessian_sign == :positive
		return 1.0
	elseif hessian_sign == :negative
		return -1.0
	elseif hessian_sign == :auto
		hessian_mw = mass_weighted_hessian(hessian, atom_masses, sign = 1.0)
		subspace = internal_subspace(atom_coords, atom_masses)

		if size(subspace, 2) == 0
			return 1.0
		end

		eigenvalues = eigvals(Symmetric(subspace' * hessian_mw * subspace))
		threshold = sqrt(eps(Float64)) * max(1.0, maximum(abs.(eigenvalues)))
		positive_count = count(>(threshold), eigenvalues)
		negative_count = count(<(-threshold), eigenvalues)

		if negative_count > positive_count
			return -1.0
		elseif positive_count > negative_count
			return 1.0
		end

		positive_weight = sum(abs, eigenvalues[eigenvalues .> threshold])
		negative_weight = sum(abs, eigenvalues[eigenvalues .< -threshold])
		return negative_weight > positive_weight ? -1.0 : 1.0
	end

	throw(ArgumentError("hessian_sign must be :auto, :positive, :negative, 1, or -1"))
end
