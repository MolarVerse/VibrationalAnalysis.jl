import VibrationalAnalysis: check_unit

"""
CLI to calculate VibrationalAnalysis on a QMCFC or a PQ input.

# Info

Calculate the wavenumbers, intensities, force constants, reduced masses and eigenvectors of a QMCFC or a PQ input.

# Args

- `restart`: The restart file.
- `hessian`: The hessian file.

# Options

- `--modesmatrix`: Write the modes in matrix notation to a file. If not specified, the modes in matrix notation will not be written.
- `--moldescriptor`: The moldescriptor file. If not specified, the intensities will not be calculated.
- `-o, --output`: The name of the output file. Default is stdout.
- `-u, --unit`: The energy unit in the Hessian file. Options are kcal, Hartree and eV. Default is kcal.

# Flags

- `--modes`: Write the modes in xyz format. If not specified, the modes will not be written.

"""
@main function vibrationalanalalysis(restart::String, hessian::String; unit = "kcal", moldescriptor = nothing, output = nothing, modes_matrix = nothing, modes::Bool = false)

	# Check if the unit is valid
	wavenumber = check_unit(unit)

	# Read restart restart file 
	atom_names, atom_masses, atom_coords, atom_types = read_rst(restart)

	# Read the hessian
	hessian = read_hessian(hessian)

	# if moldescriptor_file is not nothing
	if moldescriptor == nothing

		# calculate the wavenumbers, force constants, reduced masses and eigenvectors
		wavenumbers, force_constants, reduced_masses, eigenvectors_internal_normalized = calculate(atom_masses, atom_coords, hessian, wavenumber = wavenumber)

		# write wavenumbers, force constants, reduced masses
		write_calculate_output(wavenumbers, force_constants, reduced_masses, filename = output)
	else

		# Read the atom charges
		atom_charges = read_moldescriptor(moldescriptor, atom_names, atom_types)

		# calculate the wavenumbers, intensities, force constants, reduced masses and eigenvectors
		wavenumbers, intensities, force_constants, reduced_masses, eigenvectors_internal_normalized = calculate(atom_masses, atom_coords, atom_charges, hessian, wavenumber = wavenumber)

		write_calculate_output(wavenumbers, intensities, force_constants, reduced_masses, filename = output)
	end

	# write the modes to a file
	if modes_matrix != nothing
		write_modes(eigenvectors_internal_normalized, filename = modes_matrix)
	end

	# write modes in xyz format
	if modes
		write_modes(eigenvectors_internal_normalized, atom_coords, atom_names, filename = "modes")
	end

	return
end
