# Examples

Here are some examples of how to use the `VibrationalAnalysis` package.

Calculate the vibrational analysis of a molecule from the output of a PQ calculation.
```julia
using VibrationalAnalysis

atom_names, atom_masses, atom_coords, atom_types = read_rst("path/to/file.rst")
hessian = read_hessian("path/to/file.hessian")
atom_charges = read_moldescriptor("path/to/file.moldescriptor", atom_names, atom_types)

wavenumbers, intensities, force_constants, reduced_masses, normal_modes =
    calculate(atom_masses, atom_coords, atom_charges, hessian)

write_calculate_output(wavenumbers, intensities, force_constants, reduced_masses)
```

Calculate the vibrational analysis without intensities by providing atomic masses, atomic coordinates, and the Hessian directly.

```julia
wavenumbers, force_constants, reduced_masses, normal_modes =
    calculate(atom_masses, atom_coords, hessian)
```
