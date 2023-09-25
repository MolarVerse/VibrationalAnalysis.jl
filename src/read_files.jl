# Read Restart File

export read_rst, read_hessian, read_moldescriptor

include("massesdict.jl")

"""
    read_rst(rst_file::String)

Reads a restart file and returns a tuple of atom names, masses, and coordinates.
    
"""
function read_rst(rst_file::String)

    # Check if the file exists and 
    if !isfile(rst_file)
        error("The restart file does not exist.")
    end
    
    # Check if the file is empty
    if filesize(rst_file) == 0
        error("The restart file is empty.")
    end

    # Read the restart file
    rst_lines = readlines(rst_file)

    # Delete empty lines
    rst_lines = filter(x -> x != "", rst_lines)

    # Checks if a line begins with "Box" or "Step" and deletes these rst_lines
    # Example:
    # Step 0
    # Box  0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00
    # or 
    # Box  0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00
    # or 
    #  Step 0
    #   Box  0.0000000000000000E+00  0.0000000000000000E+00  0.0000000000000000E+00    
    rst_lines = filter(x -> !occursin(r"\s*Box", x) && !occursin(r"\s*Step", x), rst_lines)

    # Collect the atom names
    atom_names = map(x -> String(split(x)[1]), rst_lines)

    # Collect atom types
    atom_types = map(x -> parse(Int64, split(x)[3]), rst_lines)

    # Collect the masses
    atom_masses = (x -> masses[lowercase(x)]).(atom_names)

    # Collect the coordinates and convert them to Matrix{Float64}
    atom_coords = map(x -> [parse(Float64, y) for y in split(x)[4:6]], rst_lines)
    atom_coords = Matrix(hcat(atom_coords...)') # Transpose the matrix

    return atom_names, atom_masses, atom_coords, atom_types
end

"""
    read_hessian(hessian_file::String)

Reads a hessian file and returns a nxn Matrix of the hessian.

"""
function read_hessian(hessian_file::String)
    
    # Check if the file exists and
    if !isfile(hessian_file)
        error("The hessian file does not exist.")
    end

    # Check if the file is empty
    if filesize(hessian_file) == 0
        error("The hessian file is empty.")
    end

    # Delete empty lines
    hessian_lines = filter(x -> x != "", readlines(hessian_file))

    # Check if hessian_lines dont contain strings that are not numbers
    for line in hessian_lines
        for word in split(line)
            if !occursin(r"[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?", word)
                error("The hessian file contains strings that are not numbers.")
            end
        end
    end

    # Convert it to a Matrix of Float64
    hessian = [parse.(Float64, split(x)) for x in hessian_lines]
    hessian = Matrix(hcat(hessian...)') # Transpose the matrix
    
    # Check if the hessian is a square matrix
    if size(hessian)[1] != size(hessian)[2]
        error("The hessian is not a square matrix.")
    end

    return hessian

end

"""
    read_moldescriptor(moldescriptor_file::String)

Reads a moldescriptor file and returns atom_charges.
    
"""
function read_moldescriptor(moldescriptor_file::String, atom_names::Vector{String}, atom_types::Vector{Int64})

    # Check if the file exists and
    if !isfile(moldescriptor_file)
        error("The moldescriptor file does not exist.")
    end

    # Check if the file is empty
    if filesize(moldescriptor_file) == 0
        error("The moldescriptor file is empty.")
    end

    # Delete empty lines, lines containing whitespaces and comments
    moldescriptor_lines = filter(x -> x != "" && occursin(r"\w+", x) && !occursin(r"#", x), readlines(moldescriptor_file))

    # Strip the moldescriptor lines that are not lines of 3 strings
    moldescriptor_lines = filter(x -> length(split(x)) == 3, moldescriptor_lines)

    mol_types = Vector()
    i = 1

    while  i <= length(moldescriptor_lines)
        number_atoms = parse(Int64, split(moldescriptor_lines[i])[2])
        molecule = moldescriptor_lines[i+1:i+number_atoms]
        dict = Dict()
        for atom in molecule
            dict[lowercase(split(atom)[1])] = parse(Float64,split(atom)[3])
        end
        push!(mol_types, dict)
        i = i + number_atoms + 1
    end

    atom_charges = [mol_types[atom_types[i]][lowercase(atom_names[i])] for i in 1:length(atom_names)]

    # Check if the number of atoms in the moldescriptor file is the same as in the restart file
    if size(atom_charges)[1] != length(atom_names)
        error("The number of atoms in the moldescriptor file is not the same as in the restart file.")
    end

    return atom_charges
    
end