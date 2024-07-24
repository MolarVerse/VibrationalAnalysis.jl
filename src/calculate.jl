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
function calculate(atom_masses::Vector{Float64}, atom_coords::Matrix{Float64}, hessian::Matrix{Float64}; wavenumber = wavenumber_kcal)

	# Symmetrize the hessian and mass weighting
	hessian_mw = mass_weighted_hessian(hessian, atom_masses)

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
function calculate(atom_masses::Vector{Float64}, atom_coords::Matrix{Float64}, atom_charges::Vector{Float64}, hessian::Matrix{Float64}; wavenumber = wavenumber_kcal)

	# Symmetrize the hessian and mass weighting
	hessian_mw = mass_weighted_hessian(hessian, atom_masses)

	# Transform in internal coordinates
	eigenvalues, eigenvectors_internal_normalized, normalization = internal_coordinates(atom_coords, atom_masses, hessian_mw)

	# Calculate observables
	wavenumbers, omega = wavenumber(eigenvalues)
	reduced_masses = reduced_mass(normalization)
	intensities = infrared_intensity(eigenvectors_internal_normalized, atom_charges, reduced_masses)
	force_constants = force_constant(omega, reduced_masses)

	return wavenumbers, intensities, force_constants, reduced_masses, eigenvectors_internal_normalized
end
