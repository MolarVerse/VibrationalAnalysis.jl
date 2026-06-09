export read_structure, read_rst, read_xyz, read_hessian, read_moldescriptor

"""
	read_rst(rst_file::String) -> atom_names::Vector{String}, atom_masses::Vector{Float64}, atom_coords::Matrix{Float64}, atom_types::Vector{Int64}

Reads a restart file and returns a tuple of atom names, masses, and coordinates.

# Arguments
- `rst_file::String`: The restart file.
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
	read_xyz(xyz_file::String) -> atom_names::Vector{String}, atom_masses::Vector{Float64}, atom_coords::Matrix{Float64}, atom_types::Vector{Int64}

Reads a standard single-structure XYZ file and returns a tuple of atom names,
masses, coordinates, and atom types.

# Arguments
- `xyz_file::String`: The XYZ file.
"""
function read_xyz(xyz_file::String)

	# Check if the file exists and
	if !isfile(xyz_file)
		error("The xyz file does not exist.")
	end

	# Check if the file is empty
	if filesize(xyz_file) == 0
		error("The xyz file is empty.")
	end

	# Read the xyz file
	xyz_lines = readlines(xyz_file)

	# Check if the first line contains the number of atoms
	if isempty(xyz_lines)
		error("The xyz file is empty.")
	end

	number_atoms = try
		parse(Int64, strip(xyz_lines[1]))
	catch
		error("The xyz file does not start with the number of atoms.")
	end

	if number_atoms <= 0
		error("The xyz file must contain at least one atom.")
	end

	# XYZ files contain an atom-count line, a comment line and one line per atom
	if length(xyz_lines) < number_atoms + 2
		error("The xyz file does not contain the expected number of atom lines.")
	end

	atom_lines = xyz_lines[3:number_atoms+2]
	trailing_lines = xyz_lines[number_atoms+3:end]

	if any(isempty(strip(line)) for line in atom_lines)
		error("The xyz file contains empty atom lines.")
	end

	if any(!isempty(strip(line)) for line in trailing_lines)
		error("The xyz file contains more than one structure.")
	end

	atom_names = Vector{String}()
	atom_masses = Vector{Float64}()
	atom_coords = Vector{Vector{Float64}}()

	for line in atom_lines
		atom_line = split(line)

		if length(atom_line) < 4
			error("The xyz file contains atom lines that are not of the form: string float float float.")
		end

		atom_name = String(atom_line[1])
		atom_key = lowercase(atom_name)

		if !haskey(masses, atom_key)
			error("The xyz file contains an unknown atom symbol.")
		end

		coords = try
			parse.(Float64, atom_line[2:4])
		catch
			error("The xyz file contains atom lines that are not of the form: string float float float.")
		end

		push!(atom_names, atom_name)
		push!(atom_masses, masses[atom_key])
		push!(atom_coords, coords)
	end

	return atom_names, atom_masses, Matrix(hcat(atom_coords...)'), ones(Int64, number_atoms)
end

function structure_format(structure_file::String; format::Symbol = :auto)
	if format == :auto
		return lowercase(splitext(structure_file)[2]) == ".xyz" ? :xyz : :rst
	elseif format in (:rst, :xyz)
		return format
	end

	error("Unsupported structure format. Use :auto, :rst or :xyz.")
end

"""
	read_structure(structure_file::String; format::Symbol = :auto) -> atom_names::Vector{String}, atom_masses::Vector{Float64}, atom_coords::Matrix{Float64}, atom_types::Vector{Int64}

Reads either a restart file or a single-structure XYZ file.

# Arguments
- `structure_file::String`: The restart or XYZ file.
- `format::Symbol`: The structure format. Use `:auto`, `:rst` or `:xyz`.
"""
function read_structure(structure_file::String; format::Symbol = :auto)
	if structure_format(structure_file, format = format) == :xyz
		return read_xyz(structure_file)
	end

	return read_rst(structure_file)
end

function moldescriptor_types(moldescriptor_file::String)

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

	# Check if the moldescriptor lines are empty
	if length(moldescriptor_lines) == 0
		error("The moldescriptor file is not a moldescriptor file.")
	end

	# Check if the moldescriptor file contains only lines of 3 strings
	# First is a string, second is an integer and third is a float
	for line in moldescriptor_lines
		if !(occursin(r"\w+", split(line)[1]) && occursin(r"\d+", split(line)[2]) && occursin(r"[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?", split(line)[3]))
			error("The moldescriptor file contains lines that are not of the form: string integer float.")
		end
	end

	mol_types = Vector{Dict{String, Float64}}()
	i = 1

	while i <= length(moldescriptor_lines)
		number_atoms = parse(Int64, split(moldescriptor_lines[i])[2])
		molecule = moldescriptor_lines[i+1:i+number_atoms]
		dict = Dict{String, Float64}()
		for atom in molecule
			dict[lowercase(split(atom)[1])] = parse(Float64, split(atom)[3])
		end
		push!(mol_types, dict)
		i = i + number_atoms + 1
	end

	return mol_types
end

function moldescriptor_molecule_count(moldescriptor_file::String)
	return length(moldescriptor_types(moldescriptor_file))
end

"""
	read_hessian(hessian_file::String) -> hessian::Matrix{Float64}

Reads a hessian file and returns a nxn Matrix of the hessian.

# Arguments
- `hessian_file::String`: The hessian file.
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
	read_moldescriptor(moldescriptor_file::String) -> atom_charges::Vector{Float64}

Reads a moldescriptor file and returns atom_charges.

# Arguments
- `moldescriptor_file::String`: The moldescriptor file.
"""
function read_moldescriptor(moldescriptor_file::String, atom_names::Vector{String}, atom_types::Vector{Int64})

	mol_types = moldescriptor_types(moldescriptor_file)

	if any(atom_type < 1 || atom_type > length(mol_types) for atom_type in atom_types)
		error("The atom types are incompatible with the moldescriptor file.")
	end

	atom_charges = [mol_types[atom_types[i]][lowercase(atom_names[i])] for i in 1:length(atom_names)]

	# Check if the number of atoms in the moldescriptor file is the same as in the restart file
	if size(atom_charges)[1] != length(atom_names)
		error("The number of atoms in the moldescriptor file is not the same as in the restart file.")
	end

	return atom_charges

end
