export write_modes, write_wavenumber_intensity, write_calculate_output

"""
	write_modes(eigenvectors_internal_normalized, coord, atom_names; filename, amplitude, step)

Write the modes to a file in xyz format.

# Arguments
- `eigenvectors_internal_normalized::Matrix{Float64}`: The eigenvectors in internal coordinates.
- `coord::Matrix{Float64}`: The coordinates of the atoms.
- `atom_names::Vector{String}`: The names of the atoms.

# Keyword Arguments
- `filename::String`: The name of the file. Default is "modes".
- `amplitude::Float64`: The amplitude of the mode. Default is 0.25.
- `step::Float64`: The step size of the mode. Default is 0.01.

# Example
```julia
julia> using VibrationalAnalysis
julia> write_modes(eigenvectors_internal_normalized, atom_coords, atom_names)
shell> cat modes-1.xyz
```
"""
function write_modes(eigenvectors_internal_normalized::Matrix{Float64}, atom_coords::Matrix{Float64}, atom_names::Vector{String}; filename = "modes", amplitude = 0.25, step = 0.01)
	for i in 1:size(eigenvectors_internal_normalized)[2]

		# Open file
		file = open("$filename-$i.xyz", "w")

		# Number of atoms in the molecule
		n_atoms = size(atom_coords)[1]

		mode = reshape(eigenvectors_internal_normalized[:, i], 3, :)'

		for (i, alpha) in enumerate(-amplitude:step:amplitude)
			println(file, n_atoms, "\n")

			for i in 1:n_atoms
				println(file, atom_names[i], "    ",
					join(atom_coords[i, :] .+ (alpha * mode[i, :]), "   "))
			end
		end
		close(file)
	end

	return nothing
end

"""
	write_modes(eigenvectors_internal_normalized)

Write the modes to a file in matrix format.

# Arguments
- `eigenvector_internal::Matrix{Float64}`: The eigenvectors in internal coordinates.

# Keyword Arguments
- `filename::String`: The name of the file. Default is stdout.

# Example
```julia
julia> using VibrationalAnalysis
julia> write_modes(eigenvectors_internal_normalized)
shell> cat modes.dat
```
"""
function write_modes(eigenvectors_internal_normalized::Matrix{Float64}; filename = nothing)
	# Open file
	if filename != nothing
		file = open(filename, "w")
	else
		file = stdout
	end

	# Write the eigenvectors to the file in a matrix format
	for i in 1:size(eigenvectors_internal_normalized)[1]
		for j in 1:size(eigenvectors_internal_normalized)[2]
			print(file, eigenvectors_internal_normalized[i, j], " ")
		end
		println(file)
	end

	if filename != nothing
		close(file)
	end

	return nothing
end

"""
	write_calculate_output(wavenumbers, intensities; filename = stdout)

Write the wavenumbers and intensities to a file.

# Arguments
- `wavenumbers::Vector{Float64}`: The wavenumbers.
- `intensities::Vector{Float64}`: The intensities.

# Keyword Arguments
- `filename::String`: The name of the file. Default is stdout

# Example
```julia
julia> using VibrationalAnalysis
julia> write_calculate_output(wavenumbers, intensities)
```
"""
function write_calculate_output(wavenumbers, intensities; filename = nothing)

	# Open file
	if filename != nothing
		file = open(filename, "w")
	else
		file = stdout
	end

	# Write wavenumbers and intensities
	println(file, "# Wavenumbers (cm-1)  Intensities (km mol-1)")
	for i in eachindex(wavenumbers)
		# format output to 8 decimal places with scientific notation
		@printf(file, "%-8.8e\t%-8.8e\n", wavenumbers[i], intensities[i])
	end

	# Close file
	if filename != nothing
		close(file)
	end

	return nothing
end

"""
	write_calculate_output(wavenumbers, intensities, force_constants, reduced_masses; filename = stdout)

Write the wavenumbers, intensities, force constants and reduced masses to a file.

# Arguments
- `wavenumbers::Vector{Float64}`: The wavenumbers.
- `intensities::Vector{Float64}`: The intensities.
- `force_constants::Vector{Float64}`: The force constants.
- `reduced_masses::Vector{Float64}`: The reduced masses.

# Keyword Arguments
- `filename::String`: The name of the file. Default is stdout

# Example
```julia
julia> using VibrationalAnalysis
julia> write_calculate_output(wavenumbers, intensities, force_constants, reduced_masses)
```
"""
function write_calculate_output(wavenumbers, intensities, force_constants, reduced_masses; filename = nothing)

	# Open file
	if filename != nothing
		file = open(filename, "w")
	else
		file = stdout
	end

	# Write wavenumbers and intensities
	println(file, "# Wavenumbers (cm-1)  Intensities (km mol-1)  Force constants (mdyn Å-1)  Reduced masses (amu)")
	for i in eachindex(wavenumbers)
		# format output to 8 decimal places with scientific notation
		@printf(file, "%-8.8e\t%-8.8e\t%-8.8e\t%-8.8e\n", wavenumbers[i], intensities[i], force_constants[i], reduced_masses[i])
	end

	# Close file
	if filename != nothing
		close(file)
	end

	return nothing
end


"""
	write_calculate_output(wavenumbers, force_constants, reduced_masses, eigenvectors_internal_normalized; filename = stdout)

Write the wavenumbers, force constants and reduced masses to a file.

# Arguments
- `wavenumbers::Vector{Float64}`: The wavenumbers.
- `force_constants::Vector{Float64}`: The force constants.
- `reduced_masses::Vector{Float64}`: The reduced masses.

# Keyword Arguments
- `filename::String`: The name of the file. Default is stdout.

# Example
```julia
julia> using VibrationalAnalysis
julia> write_calculate_output(wavenumbers, force_constants, reduced_masses)
```
"""
function write_calculate_output(wavenumbers, force_constants, reduced_masses; filename = nothing)

	# Open file
	if filename != nothing
		file = open(filename, "w")
	else
		file = stdout
	end

	# Write wavenumbers and intensities
	println(file, "# Wavenumbers (cm-1)  Force constants (mdyn Å-1)  Reduced masses (amu)")
	for i in eachindex(wavenumbers)
		# format output to 8 decimal places with scientific notation
		@printf(file, "%-8.8e\t%-8.8e\t%-8.8e\n", wavenumbers[i], force_constants[i], reduced_masses[i])
	end

	# Close file
	if filename != nothing
		close(file)
	end

	return nothing
end
