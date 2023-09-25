
export write_modes

"""
    write_modes(eigenvector_internal::Matrix{Float64}, coord::Matrix{Float64}, atom_names::Vector{String})

Write the modes to a file.
"""

function write_modes(eigenvectors_internal_normalized::Matrix{Float64}, atom_coords::Matrix{Float64}, atom_names::Vector{String}; filename="modes", amplitude=0.25, step=0.01)
    for i in 1:size(eigenvectors_internal_normalized)[1]

        # Open file
        file = open("$filename_$i.xyz", "w")

        # Number of atoms in the molecule
        n_atoms = size(atom_coords)[1]
    
        mode = reshape(eigenvectors_internal_normalized[:, i], 3, :)'
    
        for (i,α) in enumerate(-amplitude:step:amplitude)
            println(file, n_atoms, "\n")
    
            for i in 1:n_atoms
                println(file, atom_names[i] , "    ",
                join(R[i, :] .+ (α * mode[i,:]), "   " ))
            end
        end
        close(file)
    end

    return nothing
end

"""
    write_wavenumber_frequency(wavenumbers::Vector{Float64}, intensities::Vector{Float64})

Write the wavenumbers and intensities to a file.
"""

function write_wavenumbers_intensities(wavenumbers::Vector{Float64}, intensities::Vector{Float64}; filename=stdout)
    
    # Open file
    if filename != stdout
        file = open(filename, "w")
    else
        file = filename
    end

    println(file, "Wavenumbers (cm-1)    Intensities (km mol-1)")

    for i in 1:length(wavenumbers)
        println(file, wavenumbers[i], "    ", intensities[i])
    end

    # Close file
    if filename != stdout
        close(file)
    end

    return nothing
end
