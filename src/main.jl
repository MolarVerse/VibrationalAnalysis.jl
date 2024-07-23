"""
CLI to calculate VibrationalAnalysis on a QMCFC or a PQ input.

# Args

- `rst_file`: The restart file.
- `hessian_file`: The hessian file.
- `moldescriptor_file`: The moldescriptor file. (Optional)

# Options

- `-o, --output`: The name of the output file. Default is stdout.
- `--modes`: The name of the modes file. Default is modes.dat.
- `-u, --unit`: The unit of the wavenumbers. Options are kcal, hartree and eV. Default is kcal.

"""
@main function vibanal(rst_file::String, hessian_file::String, moldescriptor_file::String; output = "stdout", modes = "modes.dat", unit = "kcal")

	# Check if the unit is valid - lowercase
	unit = lowercase(unit)
	if unit != "kcal" && unit != "hartree" && unit != "ev"
		@error "Invalid unit. Options are kcal, hartree and eV."
		return
	else
		if unit == "kcal"
			wavenumber = wavenumber_kcal
		elseif unit == "hartree"
			wavenumber = wavenumber_hartree
		elseif unit == "ev"
			wavenumber = wavenumber_eV
		else
			@error "Invalid unit. Options are kcal, hartree and eV."
			return
		end
	end

	# Read the restart file, the hessian and the atom charges
	wavenumbers, intensities, force_constants, reduced_masses, eigenvectors_internal_normalized = read_calculate(rst_file, hessian_file, moldescriptor_file, wavenumber = wavenumber)

	# Write the wavenumbers and intensities
	write_calculate_output(wavenumbers, intensities, force_constants, reduced_masses, filename = output)

	return 
end