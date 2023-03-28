using LinearAlgebra

export symmetrize_addition, symmetrize_multiplication

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
