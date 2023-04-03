
export write_modes

"""
    write_modes(eigenvector_internal::Matrix{Float64}, coord::Matrix{Float64}, atom_names::Vector{String})

Write the modes to a file.

    Parameters
    ----------
    eigenvector_internal : Matrix{Float64}
        The eigenvectors of the hessian in internal coordinates.
    coord : Matrix{Float64}
        The coordinates of the atoms.
    atom_names : Vector{String}
        The names of the atoms.
    
    Keyword Arguments
    -----------------
    filename : String
        The name of the file.

    Returns
    -------
    Nothing
"""

function write_modes(eigenvector_internal::Matrix{Float64}, coord::Matrix{Float64}, atom_names::Vector{String}; filename="modes")
    for i in 1:size(eigenvector_internal)[1]

        # Open file
        file = open("$filename_$i.xyz", "w")

        # Number of atoms in the molecule
        n_atoms = size(coord)[1]
    
        mode = reshape(eigenvector_internal[:, i], 3, :)'
    
        for (i,α) in enumerate(-0.25:0.01:0.25)
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
    write_wavenumber_frequency(wavenumber::Vector{Float64}, frequency::Vector{Float64})

Write the wavenumbers and frequencies to a file.

    Parameters
    ----------
    wavenumber : Vector{Float64}
        The wavenumbers.
    frequency : Vector{Float64}
        The frequencies.
    
    Keyword Arguments
    -----------------
    filename : String
        The name of the file.

    Returns
    -------
    Nothing
"""

function write_wavenumber_frequency(wavenumber::Vector{Float64}, frequency::Vector{Float64}; filename=stdout)
    
    # Open file
    if filename != stdout
        file = open(filename, "w")
    else
        file = filename
    end


    println(file, "Wavenumber (cm-1)    Frequency (a.u.)")

    for i in 1:length(wavenumber)
        println(file, wavenumber[i], "    ", frequency[i])
    end

    # Close file
    if filename != stdout
        close(file)
    end

    return nothing
end
