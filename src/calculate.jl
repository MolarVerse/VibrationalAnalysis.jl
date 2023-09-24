export calculate

function calculate(rst_file::String, hessian_file::String)

    # Read restart restart file 
    atom_names, atom_masses, atom_coords = read_rst(rst_file)

    # Read the hessian
    hessian = read_hessian(hessian_file)

    # Symmetrize the hessian and mass weighting
    hessian_mw = mass_weighted_hessian(hessian, atom_masses)

    # Transforme in internal coordinates
    eigenvalues, eigenvectors_internal_normalized, N = internal_coordinates(atom_coords, atom_masses, hessian_mw)

    #TODO: Build read moldescriptor or charges reader
    atom_charges = [-2.0, 1.0, 1.0]

    # Calculate observables
    wavenumbers = wavenumber_kcal(eigenvalues)
    intensities = infrared_intensity(N, atom_coords, atom_charges)
    reduced_mass = reduced_mass(N)
    force_constant = force_constant(wavenumbers, reduced_mass)
    
    # Write the modes as xyz-files
    write_modes(eigenvectors_internal_normalized, atom_coords, atom_names)
    # Write the wavenumbers and intensities to a file
    write_wavenumber_frequency(wavenumbers, intensities)

end
