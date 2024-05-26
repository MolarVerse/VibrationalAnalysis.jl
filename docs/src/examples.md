# Examples

Here are some examples of how to use the `VibrationalAnalysis` package.

Calculate the vibrational analysis of a molecule from the output of a PQ calculation.
```julia
read_calculate("path/to/file.rst", "path/to/file.hessian", "path/to/file.moldescriptor")
```

Calculate the vibrational analysis by providing the necessary data. Atomic masses, atomic coordinates, atomic charges, and the hessian are the necessary data to calculate the vibrational analysis.

```julia
calculate(atom_masses, atom_coords, atom_charges, hessian)
```