
export translational_modes, rotational_modes, transformation_matrix

"""
    translational_modes(masses::Vector{Float64})
    
    Calculate the translational modes of a molecule.
    
    Parameters
    ----------
    masses : Vector{Float64}
        The masses of the atoms.
    
    Returns
    -------
    translation : Matrix{Float64}
        The translational modes of the molecule.
    
"""

function translational_modes(masses::Vector{Float64})
    # Create translation matrix 3N x 3
    masses_diagonal = [diagm(0 => [i, i, i]) for i in masses]
    
    # Convert to Matrix{Float64}
    translation = hcat(masses_diagonal...)'

    # Normalize the translation matrix
    translation = translation ./ sqrt.(sum(translation.^2, dims=1))

    return translation
end

"""
    rotational_modes(coord::Matrix{Float64}, masses::Vector{Float64})
    
    Calculate the rotational modes of a molecule.
    
    Parameters
    ----------
    coord : Matrix{Float64}
        The coordinates of the atoms.
    masses : Vector{Float64}
        The masses of the atoms.
    
    Returns
    -------
    rotation : Matrix{Float64}
        The rotational modes of the molecule.
    
"""

function rotational_modes(coord::Matrix{Float64}, masses::Vector{Float64})
    # Translate the coordinates to the center of mass
    coord = center_to_com(coord, masses)
    # Calculate the inertia tensor
    inertia_tensor = calculate_inertia_tensor(coord, masses)
    # Calculate the eigenvalues and eigenvectors of the inertia tensor
    _ , eigenvectors = eigen(inertia_tensor)
    
    # Notation from Gaussian Vibrational Analysis
    X = eigenvectors
    # Calculate the rotational frame
    P = coord * eigenvectors
    
    # Initialize the rotation matrix 3N x 3
    rotation = zeros(size(coord, 1) * 3, 3)

    # Calculate the rotation matrix
    rotation[:, 1] += ((P[:,2] * X[3, :]' .- P[:,3] * X[2,:]') .* sqrt.(masses))[:]
    rotation[:, 2] += ((P[:,3] * X[1, :]' .- P[:,1] * X[3,:]') .* sqrt.(masses))[:]
    rotation[:, 3] += ((P[:,1] * X[2, :]' .- P[:,2] * X[1,:]') .* sqrt.(masses))[:]
    
    # Normalize the rotation matrix
    rotation = rotation ./ sqrt.(sum(rotation.^2, dims=1))
    return rotation
end

"""
    transformation_matrix(coord::Matrix{Float64}, masses::Vector{Float64})
    
    Calculate the transformation matrix.
    
    Parameters
    ----------
    coord : Matrix{Float64}
        The coordinates of the atoms.
    masses : Vector{Float64}
        The masses of the atoms.
    
    Returns
    -------
    transformation : Matrix{Float64}
        The transformation matrix.
    
"""

function transformation_matrix(coord::Matrix{Float64}, masses::Vector{Float64})
    
    # Initialize transformation matrix 3N x 3N
    transformation = zeros(size(coord, 1) * 3, size(coord, 1) * 3)
    
    # Calculate the translational modes
    translation = translational_modes(masses)
    # Calculate the rotational modes
    rotation = rotational_modes(coord, masses)
    
    # Combine the translational and rotational modes
    transformation[:, 1:3] = translation
    transformation[:, 4:6] = rotation

    return transformation
end