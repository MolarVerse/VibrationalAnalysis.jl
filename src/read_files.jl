# Read Restart File

export read_rst, read_hessian

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
    atom_names = map(x -> split(x)[1], rst_lines)

    # Collect the masses
    atom_masses = (x->masses[lowercase(x)]).(atom_names)

    # Collect the coordinates and convert them to Matrix{Float64}
    atom_coords = map(x -> [parse(Float64, y) for y in split(x)[4:6]], rst_lines)
    atom_coords = Matrix(hcat(atom_coords...)') # Transpose the matrix

    return atom_names, atom_masses, atom_coords
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

    # Convert it to a Matrix of Float64
    hessian = [parse.(Float64, split(x)) for x in hessian_lines]
    hessian = Matrix(hcat(hessian...)') # Transpose the matrix
    

    # Check if the hessian is a square matrix
    if size(hessian)[1] != size(hessian)[2]
        error("The hessian is not a square matrix.")
    end

    return hessian

end