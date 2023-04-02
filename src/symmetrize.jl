using LinearAlgebra

export symmetrize_addition, symmetrize_multiplication, mass_weighted_hessian

"""
    symmetrize_addition(hessian::Matrix{Float64})

Symmetrize a Hessian matrix by adding the transpose and dividing by 2.
"""

function symmetrize_addition(hessian::Matrix{Float64})
    return (hessian + hessian') / 2
end

"""
    symmetrize_multiplication(hessian::Matrix{Float64})

Symmetrize a Hessian matrix by multiplying by the transpose and taking the square root of the absolute value of the eigenvalues.
"""

function symmetrize_multiplication(hessian::Matrix{Float64})
    _H = hessian' * hessian
    v, V = eigen(_H)
    _H_sym = V * Diagonal(sqrt.(abs.(v))) * V'
    return _H_sym
end

"""
    mass_weight_hessian(hessian::Matrix{Float64}, atom_masses::Vector{Float64})

Mass weight a Hessian matrix.
"""

function mass_weighted_hessian(hessian::String)
    # Read the hessian
    hessian = read_hessian(hessian)
    # Read the rst file
    hessian = symmetrize_multiplication(hessian)
    
    # Create 3N x 3N matrix of masses
    masses_repeat = repeat(masses, inner=3)
    # Matrix of masses \sqrt{m_i m_j}
    masses_matrix = (masses_repeat * masses_repeat').^(1/2)

    # Mass weight the hessian
    hessian = hessian ./ masses_mat

    return hessian
end