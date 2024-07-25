"""
	symmetrize_addition(hessian)

Symmetrize a Hessian matrix by adding the transpose and dividing by 2.

``H_sym = (H + H') / 2``

# Arguments
- `hessian::Matrix{Float64}`: The Hessian matrix.

# Returns
- `hessian_sym::Matrix{Float64}`: The symmetrized Hessian matrix.
"""
function symmetrize_addition(hessian::Matrix{Float64})
	return (hessian + hessian') / 2
end

"""
	symmetrize_multiplication(hessian)

Symmetrize a Hessian matrix by multiplying by the transpose and taking the square root of the absolute value of the eigenvalues.

``H_2 = H * H'``

``v, V = eigen(H_2)``

``H_sym = V * √|v| * V'``

# Arguments
- `hessian::Matrix{Float64}`: The Hessian matrix.

# Returns
- `hessian_sym::Matrix{Float64}`: The symmetrized Hessian matrix.
"""
function symmetrize_multiplication(hessian::Matrix{Float64})
	_H = hessian' * hessian
	v, V = eigen(_H)
	_H_sym = V * Diagonal(sqrt.(abs.(v))) * V'
	return _H_sym
end

"""
	mass_weight_hessian(hessian, atom_masses)

Mass weight a Hessian matrix using matrix multiplication for symmetrization.

# Arguments
- `hessian::Matrix{Float64}`: The Hessian matrix.
- `atom_masses::Vector{Float64}`: The masses of the atoms.

# Returns
- `hessian::Matrix{Float64}`: The mass weighted Hessian matrix.
"""
function mass_weighted_hessian(hessian::Matrix{Float64}, atom_masses::Vector{Float64})
	# Symmetrize the hessian
	hessian = symmetrize_multiplication(hessian)

	# Get masses matrix
	_masses_matrix = masses_matrix(atom_masses)

	# Mass weight the hessian
	hessian = hessian ./ _masses_matrix

	return hessian
end

"""
	mass_weight_hessian_add(hessian, atom_masses)

Mass weight a Hessian matrix using symmetrize by addition.

# Arguments
- `hessian::Matrix{Float64}`: The Hessian matrix.
- `atom_masses::Vector{Float64}`: The masses of the atoms.

# Returns
- `hessian::Matrix{Float64}`: The mass weighted Hessian matrix.
"""
function mass_weighted_hessian_add(hessian::Matrix{Float64}, atom_masses::Vector{Float64})
	# Symmetrize the hessian
	hessian = symmetrize_addition(hessian)

	# Get masses matrix
	_masses_matrix = masses_matrix(atom_masses)

	# Mass weight the hessian
	hessian = hessian ./ _masses_matrix

	return hessian
end

"""
	masses_matrix(atom_masses)

Create a matrix of masses from a vector of atom masses.

# Arguments
- `atom_masses::Vector{Float64}`: The masses of the atoms.

# Returns
- `masses_matrix::Matrix{Float64}`: The matrix of masses.
"""
function masses_matrix(atom_masses::Vector{Float64})
	
	# Create 3N x 3N matrix of masses
	masses_repeat = repeat(atom_masses, inner = 3)
	# Matrix of masses \sqrt{m_i m_j}
	_masses_matrix = (masses_repeat * masses_repeat') .^ (1 / 2)

	return _masses_matrix
end
