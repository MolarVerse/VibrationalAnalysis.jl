# Read Restart File

export read_rst

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
    
    # Read the restart file
    rst_lines = readlines(rst_file)

    # Checks if the file is empty
    if rst_lines == ""
        error("The restart file is empty.")
    end

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
    atom_coords = hcat(atom_coords...)

    return atom_names, atom_masses, atom_coords

end