
export center_to_com, inertia_tensor

"""
    center_to_com(coord::Matrix{Float64}, masses::Vector{Float64})

Translate the coordinates to the center of mass.
    
    Parameters
    ----------
    coord : Matrix{Float64}
        The coordinates of the atoms.
    masses : Vector{Float64}
        The masses of the atoms.
    
    Returns
    -------
    coord : Matrix{Float64}
        The coordinates of the atoms translated to the center of mass.

"""

function center_to_com(coord::Matrix{Float64}, masses::Vector{Float64})
    # Calculate the center of mass
    com = sum(coord .* masses, dims=1) / sum(masses)
    # Translate the coordinates to the center of mass
    coord = coord .- com
    return coord
end

"""
    inertia_tensor(coord::Matrix{Float64}, masses::Vector{Float64})

Calculate the inertia tensor.
    
        Parameters
        ----------
        coord : Matrix{Float64}
            The coordinates of the atoms.
        masses : Vector{Float64}
            The masses of the atoms.
    
        Returns
        -------
        inertia_tensor : Matrix{Float64}
            The inertia tensor.
    
    """

function inertia_tensor(coord::Matrix{Float64}, masses::Vector{Float64})
    
    # center the coordinates to the center of mass
    coord = center_to_com(coord, masses)

    # Calculate the inertia tensor
    x = coord[:, 1]
    y = coord[:, 2]
    z = coord[:, 3]
    Ixx = sum(masses .* (y.^2 + z.^2))
    Iyy = sum(masses .* (x.^2 + z.^2))
    Izz = sum(masses .* (x.^2 + y.^2))
    Ixy = - sum(masses .* (x .* y))
    Ixz = - sum(masses .* (x .* z))
    Iyz = - sum(masses .* (y .* z))

    return [Ixx Ixy Ixz; Ixy Iyy Iyz; Ixz Iyz Izz]
end
