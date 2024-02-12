"""
    center_to_com(coord::Matrix{Float64}, masses::Vector{Float64})

Translate the coordinates to the center of mass.

# Arguments
- `coord::Matrix{Float64}`: The coordinates.
- `masses::Vector{Float64}`: The masses.

# Example
```julia-repl
julia> center_to_com(coord, masses)
```
"""
function center_to_com(atom_coords::Matrix{Float64}, atom_masses::Vector{Float64})
    
    # Calculate the center of mass
    com = sum(atom_coords .* atom_masses, dims=1) / sum(atom_masses)
    
    # Translate the coordinates to the center of mass
    atom_coords = atom_coords .- com#

    return atom_coords
end

"""
    inertia_tensor(coord::Matrix{Float64}, masses::Vector{Float64})

Calculate the inertia tensor of a molecule.

# Arguments
- `coord::Matrix{Float64}`: The coordinates.
- `masses::Vector{Float64}`: The masses.

# Example
```julia-repl
julia> inertia_tensor(coord, masses)
```
"""
function inertia_tensor(atom_coords::Matrix{Float64}, atom_masses::Vector{Float64})
    
    # center the coordinates to the center of mass
    atom_coords = center_to_com(atom_coords, atom_masses)

    # Calculate the inertia tensor
    x = atom_coords[:, 1]
    y = atom_coords[:, 2]
    z = atom_coords[:, 3]
    Ixx = sum(atom_masses .* (y.^2 + z.^2))
    Iyy = sum(atom_masses .* (x.^2 + z.^2))
    Izz = sum(atom_masses .* (x.^2 + y.^2))
    Ixy = - sum(atom_masses .* (x .* y))
    Ixz = - sum(atom_masses .* (x .* z))
    Iyz = - sum(atom_masses .* (y .* z))

    return [Ixx Ixy Ixz; Ixy Iyy Iyz; Ixz Iyz Izz]
end
