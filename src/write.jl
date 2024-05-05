export write_modes, write_wavenumber_intensity

"""
    write_modes(eigenvector_internal::Matrix{Float64}, coord::Matrix{Float64}, atom_names::Vector{String}; filename="modes", amplitude=0.25, step=0.01) -> nothing

Write the modes to a file in xyz format.

# Arguments
- `eigenvector_internal::Matrix{Float64}`: The eigenvectors in internal coordinates.
- `coord::Matrix{Float64}`: The coordinates of the atoms.
- `atom_names::Vector{String}`: The names of the atoms.

# Keyword Arguments
- `filename::String`: The name of the file. Default is "modes".
- `amplitude::Float64`: The amplitude of the mode. Default is 0.25.
- `step::Float64`: The step size of the mode. Default is 0.01.
"""
function write_modes(eigenvectors_internal_normalized::Matrix{Float64}, atom_coords::Matrix{Float64}, atom_names::Vector{String}; filename="modes", amplitude=0.25, step=0.01)
    for i in 1:size(eigenvectors_internal_normalized)[2]

        # Open file
        file = open("$filename-$i.xyz", "w")

        # Number of atoms in the molecule
        n_atoms = size(atom_coords)[1]
    
        mode = reshape(eigenvectors_internal_normalized[:, i], 3, :)'
    
        for (i,alpha) in enumerate(-amplitude:step:amplitude)
            println(file, n_atoms, "\n")
    
            for i in 1:n_atoms
                println(file, atom_names[i] , "    ",
                join(atom_coords[i, :] .+ (alpha * mode[i,:]), "   " ))
            end
        end
        close(file)
    end

    return nothing
end

"""
    write_wavenumber_frequency(wavenumbers::Vector{Float64}, intensities::Vector{Float64}; filename=stdout) -> nothing

Write the wavenumbers and intensities to a file.

# Arguments
- `wavenumbers::Vector{Float64}`: The wavenumbers.
- `intensities::Vector{Float64}`: The intensities.

# Keyword Arguments
- `filename::String`: The name of the file. Default is stdout.
"""
function write_wavenumber_intensity(wavenumbers, intensities; filename=stdout)
    
    # Open file
    if filename != stdout
        file = open(filename, "w")
    else
        file = filename
    end

    println(file, "# Wavenumbers (cm-1)    Intensities (km mol-1)")

    for i in eachindex(wavenumbers)
        println(file, wavenumbers[i], "    ", intensities[i])
    end

    # Close file
    if filename != stdout
        close(file)
    end

    return nothing
end
