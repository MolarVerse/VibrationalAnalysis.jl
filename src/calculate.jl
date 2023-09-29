export calculate

function calculate(rst_file::String, hessian_file::String, moldescriptor_file::String)

    # Read restart restart file 
    atom_names, atom_masses, atom_coords, atom_types = read_rst(rst_file)

    # Read the hessian
    hessian = read_hessian(hessian_file)

    # Symmetrize the hessian and mass weighting
    hessian_mw = mass_weighted_hessian(hessian, atom_masses)

    # Transforme in internal coordinates
    eigenvalues, eigenvectors_internal_normalized, normalization = internal_coordinates(atom_coords, atom_masses, hessian_mw)

    #Read the atom charges
    atom_charges = read_moldescriptor(moldescriptor_file, atom_names, atom_types)

    # Calculate observables
    wavenumbers = wavenumber_kcal(eigenvalues)
    reduced_masses = reduced_mass(normalization)
    intensities = infrared_intensity(eigenvectors_internal_normalized, atom_charges, reduced_masses)
    force_constants = force_constant(wavenumbers, reduced_masses)
    
    # Write the modes as xyz-files
    write_modes(eigenvectors_internal_normalized, atom_coords, atom_names)
    # Write the wavenumbers and intensities to a file
    write_wavenumber_intensity(wavenumbers, intensities)

end
