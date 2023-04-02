
export center_to_com, calculate_inertia_tensor

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
    calculate_inertia_tensor(coord::Matrix{Float64}, masses::Vector{Float64})

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

function calculate_inertia_tensor(coord::Matrix{Float64}, masses::Vector{Float64})
    
    # Initialize the inertia tensor
    inertia_tensor = zeros(3, 3)

    # Calculate the inertia tensor
    for i in 1:size(coord, 1)
        inertia_tensor += masses[i] .* (coord[i, :] .* coord[i, :]' .- coord[i, :] .* coord[i, :])
    end

    return inertia_tensor
end
