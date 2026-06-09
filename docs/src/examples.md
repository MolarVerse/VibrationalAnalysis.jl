# Examples

Here are some examples of how to use the `VibrationalAnalysis` package.

Calculate the vibrational analysis of a molecule from a structure file.
```julia
using VibrationalAnalysis

atom_names, atom_masses, atom_coords, atom_types = read_structure("path/to/file.rst")
hessian = read_hessian("path/to/file.hessian")
atom_charges = read_moldescriptor("path/to/file.moldescriptor", atom_names, atom_types)

wavenumbers, intensities, force_constants, reduced_masses, normal_modes =
    calculate(atom_masses, atom_coords, atom_charges, hessian)

write_calculate_output(wavenumbers, intensities, force_constants, reduced_masses)
```

The same workflow works with a standard single-structure XYZ file. If intensities are needed, the moldescriptor file must contain exactly one molecule definition.

```julia
using VibrationalAnalysis

atom_names, atom_masses, atom_coords, atom_types = read_structure("path/to/file.xyz")
hessian = read_hessian("path/to/file.hessian")

wavenumbers, force_constants, reduced_masses, normal_modes =
    calculate(atom_masses, atom_coords, hessian)
```

Calculate the vibrational analysis without intensities by providing atomic masses, atomic coordinates, and the Hessian directly.

```julia
wavenumbers, force_constants, reduced_masses, normal_modes =
    calculate(atom_masses, atom_coords, hessian)
```
