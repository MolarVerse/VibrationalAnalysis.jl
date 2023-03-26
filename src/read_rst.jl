# Read Restart File

export read_rst

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

    # Checks if the file has "Box" included in the first, second or last LinearAlgebra
    if rst_lines[1:3] == "Box" || rst_lines[4:6] == "Box" || rst_lines[end-2:end] == "Box"
        error("The restart file is not in the correct format.")
    end
    
end