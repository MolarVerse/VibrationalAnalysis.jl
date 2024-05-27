"""
	symmetrize_addition(hessian::Matrix{Float64}) -> hessian_sym::Matrix{Float64}

Symmetrize a Hessian matrix by adding the transpose and dividing by 2.

``H_sym = (H + H') / 2``

# Arguments
- `hessian::Matrix{Float64}`: The Hessian matrix.
"""
function symmetrize_addition(hessian::Matrix{Float64})
	return (hessian + hessian') / 2
end

"""
	symmetrize_multiplication(hessian::Matrix{Float64}) -> hessian_sym::Matrix{Float64}

Symmetrize a Hessian matrix by multiplying by the transpose and taking the square root of the absolute value of the eigenvalues.

``H_2 = H * H'``

``v, V = eigen(H_2)``

``H_sym = V * √|v| * V'``

# Arguments
- `hessian::Matrix{Float64}`: The Hessian matrix.
"""
function symmetrize_multiplication(hessian::Matrix{Float64})
	_H = hessian' * hessian
	v, V = eigen(_H)
	_H_sym = V * Diagonal(sqrt.(abs.(v))) * V'
	return _H_sym
end

"""
	mass_weight_hessian(hessian::Matrix{Float64}, atom_masses::Vector{Float64}) -> hessian::Matrix{Float64}

Mass weight a Hessian matrix using matrix multiplication for symmetrization.

# Arguments
- `hessian::Matrix{Float64}`: The Hessian matrix.
- `atom_masses::Vector{Float64}`: The masses of the atoms.
"""
function mass_weighted_hessian(hessian::Matrix{Float64}, atom_masses::Vector{Float64})
	# Symmetrize the hessian
	hessian = symmetrize_multiplication(hessian)

	# Create 3N x 3N matrix of masses
	masses_repeat = repeat(atom_masses, inner = 3)
	# Matrix of masses \sqrt{m_i m_j}
	masses_matrix = (masses_repeat * masses_repeat') .^ (1 / 2)

	# Mass weight the hessian
	hessian = hessian ./ masses_matrix

	return hessian
end

"""
	mass_weight_hessian_add(hessian::Matrix{Float64}, atom_masses::Vector{Float64})

Mass weight a Hessian matrix using symmetrize by addition.

# Arguments
- `hessian::Matrix{Float64}`: The Hessian matrix.
- `atom_masses::Vector{Float64}`: The masses of the atoms.
"""
function mass_weighted_hessian_add(hessian::Matrix{Float64}, atom_masses::Vector{Float64})
	# Symmetrize the hessian
	hessian = symmetrize_addition(hessian)

	# Create 3N x 3N matrix of masses
	masses_repeat = repeat(atom_masses, inner = 3)
	# Matrix of masses \sqrt{m_i m_j}
	masses_matrix = (masses_repeat * masses_repeat') .^ (1 / 2)

	# Mass weight the hessian
	hessian = hessian ./ masses_matrix

	return hessian
end
