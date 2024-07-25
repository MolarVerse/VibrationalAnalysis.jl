"""
	translational_modes(atom_masses)
	
Calculate the translational modes of a molecule. 

# Arguments
- `atom_masses::Vector{Float64}`: The masses of the atoms.

# Returns
- `translation::Matrix{Float64}`: The translational modes.
"""
function translational_modes(atom_masses::Vector{Float64})
	# Create translation matrix 3N x 3
	masses_diagonal = [diagm(0 => [i, i, i]) for i in atom_masses]

	# Convert to Matrix{Float64}
	translation = sqrt.(hcat(masses_diagonal...)')

	# Normalize the translation matrix
	translation = translation ./ sqrt.(sum(translation .^ 2, dims = 1))

	return translation
end

"""
	rotational_modes(atom_coords, atom_masses)
	
Calculate the rotational modes of a molecule.

# Arguments
- `atom_coords::Matrix{Float64}`: The coordinates of the atoms.
- `atom_masses::Vector{Float64}`: The masses of the atoms.

# Returns
- `rotation::Matrix{Float64}`: The rotational modes.
"""
function rotational_modes(atom_coords::Matrix{Float64}, atom_masses::Vector{Float64})

	# Translate the coordinates to the center of mass
	atom_coords_cm = center_to_com(atom_coords, atom_masses)
	# Calculate the inertia tensor
	inertia_tens = inertia_tensor(atom_coords, atom_masses)
	# Calculate the eigenvalues and eigenvectors of the inertia tensor
	_, eigenvectors = eigen(inertia_tens)

	# Notation from Gaussian Vibrational Analysis
	X = eigenvectors
	# Calculate the rotational frame
	P = atom_coords_cm * eigenvectors

	# Initialize the rotation matrix 3N x 3
	rotation = zeros(size(atom_coords, 1) * 3, 3)

	# Calculate the rotation matrix
	rotation[:, 1] += ((P[:, 2]*X[3, :]'.-P[:, 3]*X[2, :]').*sqrt.(atom_masses))[:]
	rotation[:, 2] += ((P[:, 3]*X[1, :]'.-P[:, 1]*X[3, :]').*sqrt.(atom_masses))[:]
	rotation[:, 3] += ((P[:, 1]*X[2, :]'.-P[:, 2]*X[1, :]').*sqrt.(atom_masses))[:]

	# Normalize the rotation matrix
	rotation = rotation ./ sqrt.(sum(rotation .^ 2, dims = 1))

	return rotation
end

"""
	transformation_matrix(atom_coords, atom_masses)
	
Calculate the transformation matrix.

# Arguments
- `atom_coords::Matrix{Float64}`: The coordinates of the atoms.
- `atom_masses::Vector{Float64}`: The masses of the atoms.

# Returns
- `transformation::Matrix{Float64}`: The transformation matrix.
"""
function transformation_matrix(atom_coords::Matrix{Float64}, atom_masses::Vector{Float64})

	# Initialize transformation matrix 3N x 6
	transformation = zeros(size(atom_coords, 1) * 3, 6)

	# Calculate the translational modes
	translation = translational_modes(atom_masses)
	# Calculate the rotational modes
	rotation = rotational_modes(atom_coords, atom_masses)

	# Combine the translational and rotational modes
	transformation[:, 1:3] = translation
	transformation[:, 4:6] = rotation

	return transformation
end

"""
	internal_coordinates(atom_coords, atom_masses, hessian_mw)
	
Calculate the internal coordinates of a molecule.

# Arguments
- `atom_coords::Matrix{Float64}`: The coordinates of the atoms.
- `atom_masses::Vector{Float64}`: The masses of the atoms.
- `hessian_mw::Matrix{Float64}`: The mass weighted hessian.

# Returns
- `eigenvalues::Vector{Float64}`: The eigenvalues of the hessian in internal coordinates.
- `eigenvectors_internal_normalized::Matrix{Float64}`: The eigenvectors in internal coordinates.
- `normalization::Vector{Float64}`: The normalization of the eigenvectors.
"""
function internal_coordinates(atom_coords::Matrix{Float64}, atom_masses::Vector{Float64}, hessian_mw::Matrix{Float64})

	# Calculate the transformation matrix
	transformation = transformation_matrix(atom_coords, atom_masses)
	# Gram-Schmidt orthogonalization
	transformation = qr(transformation).Q
	# Calculate the hessian in internal coordinates
	internal_hessian = transformation' * hessian_mw * transformation
	# Calculate the eigenvalues and eigenvectors of the hessian in internal coordinates
	eigenvalues, eigenvectors = eigen(internal_hessian)

	# Calculate the masses matrix
	masses_matrix = diagm(0 => [1 / sqrt(i) for i in atom_masses for j in 1:3])
	# Calculate the eigenvectors in internal coordinates
	eigenvectors_internal = masses_matrix * transformation * eigenvectors

	# Normalize Vectors
	normalization = sqrt.(1 ./ sum(eigenvectors_internal .^ 2, dims = 1))[:] # convert to vector
	# Normalize the eigenvectors
	eigenvectors_internal_normalized = eigenvectors_internal .* normalization'

	return eigenvalues, eigenvectors_internal_normalized, normalization
end