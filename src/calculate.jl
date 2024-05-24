export calculate, read_calculate

"""
	read_calculate(rst_file::String, hessian_file::String, moldescriptor_file::String) -> wavenumbers::Vector{Float64}, intensities::Vector{Float64}, force_constants::Vector{Float64}, reduced_masses::Vector{Float64}

Reads the restart file, the hessian and the atom charges and calculates the wavenumbers, intensities, force constants and reduced masses.

# Arguments
- `rst_file::String`: The restart file.
- `hessian_file::String`: The hessian file.
- `moldescriptor_file::String`: The moldescriptor file.

# Example
```julia-repl
julia> read_calculate("restart.rst", "hessian.dat", "moldescriptor.dat")
```
"""
function read_calculate(rst_file::String, hessian_file::String, moldescriptor_file::String)

	# Read restart restart file 
	atom_names, atom_masses, atom_coords, atom_types = read_rst(rst_file)

	# Read the hessian
	hessian = read_hessian(hessian_file)

	# Read the atom charges
	atom_charges = read_moldescriptor(moldescriptor_file, atom_names, atom_types)

	wavenumbers, intensities, force_constants, reduced_masses = calculate(atom_masses, atom_coords, atom_charges, hessian)

	return wavenumbers, intensities, force_constants, reduced_masses
end

"""
	calculate(atom_masses::Vector{Float64}, atom_coords::Matrix{Float64}, atom_charges::Vector{Float64}, hessian::Matrix{Float64}) -> wavenumbers::Vector{Float64}, intensities::Vector{Float64}, force_constants::Vector{Float64}, reduced_masses::Vector{Float64}

Calculates the wavenumbers, intensities, force constants and reduced masses from the atom masses, atom coordinates, atom charges and the hessian.

# Arguments
- `atom_masses::Vector{Float64}`: Vector of atom_masses (n)
- `atom_coords::Matrix{Float64}`: Matrix of atom coordinates (3xn)
- `atom_charges::Vector{Float64}`: Vector of atom charges (n)
- `hessian::Matrix{Float64}`: Matrix of the hessian (3nx3n)

# Optional arguments
- `wavenumber::Function`: The wavenumber function to use. Either `wavenumber_kcal` or `wavenumber_dftb`. Default is `wavenumber_kcal`. 

# Example
```julia-repl
julia> calculate(atom_masses, atom_coords, atom_charges, hessian)
julia> calculate(atom_masses, atom_coords, atom_charges, hessian, wavenumber=wavenumber_dftb)
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

	return wavenumbers, intensities, force_constants, reduced_masses
end
