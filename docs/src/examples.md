# Examples

This page focuses on API-level examples beyond the quickstart.

## Restart input with intensities

```julia
using VibrationalAnalysis

atom_names, atom_masses, atom_coords, atom_types = read_structure("path/to/file.rst")
hessian = read_hessian("path/to/file.hessian")
atom_charges = read_moldescriptor("path/to/file.moldescriptor", atom_names, atom_types)

wavenumbers, intensities, force_constants, reduced_masses, normal_modes =
    calculate(atom_masses, atom_coords, atom_charges, hessian)

write_calculate_output(wavenumbers, intensities, force_constants, reduced_masses)
```

## XYZ input without intensities

The same workflow works with a standard single-structure XYZ file. If intensities are needed, the moldescriptor file must contain exactly one molecule definition.

```julia
using VibrationalAnalysis

atom_names, atom_masses, atom_coords, atom_types = read_structure("path/to/file.xyz")
hessian = read_hessian("path/to/file.hessian")

wavenumbers, force_constants, reduced_masses, normal_modes =
    calculate(atom_masses, atom_coords, hessian)
```

## Choosing the Hessian unit convention

```julia
using VibrationalAnalysis

wavenumbers_kcal, force_constants, reduced_masses, normal_modes =
    calculate(atom_masses, atom_coords, hessian)

wavenumbers_hartree, force_constants, reduced_masses, normal_modes =
    calculate(atom_masses, atom_coords, hessian, wavenumber = wavenumber_hartree)

wavenumbers_ev, force_constants, reduced_masses, normal_modes =
    calculate(atom_masses, atom_coords, hessian, wavenumber = wavenumber_eV)
```

## Controlling the Hessian sign convention

`calculate(...)` defaults to `hessian_sign = :auto`. Use an explicit choice when the source convention is known:

```julia
using VibrationalAnalysis

wavenumbers, force_constants, reduced_masses, normal_modes =
    calculate(atom_masses, atom_coords, hessian, hessian_sign = :positive)

wavenumbers, force_constants, reduced_masses, normal_modes =
    calculate(atom_masses, atom_coords, hessian, hessian_sign = :negative)

wavenumbers, force_constants, reduced_masses, normal_modes =
    calculate(atom_masses, atom_coords, hessian, hessian_sign = -1)
```

## Direct array-only workflow

You can bypass file readers entirely if your code already has masses, coordinates, charges, and Hessians in memory.

```julia
wavenumbers, force_constants, reduced_masses, normal_modes =
    calculate(atom_masses, atom_coords, hessian)
```

With intensities:

```julia
wavenumbers, intensities, force_constants, reduced_masses, normal_modes =
    calculate(atom_masses, atom_coords, atom_charges, hessian)
```

## Writing normal modes

Write the mode matrix:

```julia
write_modes(normal_modes, filename = "normal_modes.dat")
```

Write animated XYZ trajectories:

```julia
write_modes(normal_modes, atom_coords, atom_names, filename = "modes")
```
