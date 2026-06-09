# Input Files

`VibrationalAnalysis.jl` currently reads four input file types:

- restart files through `read_rst(...)`
- standard single-structure XYZ files through `read_xyz(...)`
- Hessian matrices through `read_hessian(...)`
- moldescriptor files through `read_moldescriptor(...)`

The recommended public entrypoint for structure data is `read_structure(...)`.

## Structure dispatch

```julia
atom_names, atom_masses, atom_coords, atom_types = read_structure("structure.rst")
atom_names, atom_masses, atom_coords, atom_types = read_structure("structure.xyz")
atom_names, atom_masses, atom_coords, atom_types = read_structure("structure.xyz", format = :xyz)
```

Automatic dispatch is extension-based:

- `.xyz` -> XYZ reader
- everything else -> restart reader

Use the `format` keyword when you want to make the choice explicit.

## Restart files

The restart reader expects one atom per line. In the current parser:

- column `1` is the atom symbol
- column `3` is the atom type
- columns `4:6` are Cartesian coordinates
- blank lines are ignored
- lines starting with `Step` or `Box` are ignored

Example:

```text
Step 100
Box 10000.0 10000.0 10000.0
o   1   1   -0.01294293652  -0.00022738404  -0.00072569662
h   1   1    0.43396012267   0.59193103593   0.66981595463
h   1   1    0.68120631717  -0.49583217141  -0.52277201300
```

Restart files carry `atom_types`, so they can be paired with multi-molecule moldescriptor files.

## XYZ files

The XYZ reader expects standard single-structure XYZ:

- first line: number of atoms
- second line: comment
- remaining lines: `symbol x y z`

Example:

```text
3
water
O -0.01294293652 -0.00022738404 -0.00072569662
H 0.43396012267 0.59193103593 0.66981595463
H 0.68120631717 -0.49583217141 -0.52277201300
```

Current constraints:

- multi-frame XYZ files are rejected
- unknown atom symbols are rejected
- `atom_types` default to `1` for every atom

That last point matters for intensities: XYZ files do not encode molecule-type assignments, so only single-molecule moldescriptor files are compatible with XYZ input.

## Moldescriptor files

`read_moldescriptor(...)` maps atom symbols and atom types to partial charges.

The parser looks for repeated molecule blocks of the form:

```text
H2O 3 0.0
O 0 -0.65966
H 1 0.32983
H 1 0.32983
```

The current implementation:

- ignores comments and blank lines outside the data blocks
- reads one charge per atom
- validates that each `atom_type` points to an available molecule block

With restart input, `atom_types` can select among multiple molecule definitions. With XYZ input, all atoms have type `1`, so the moldescriptor file must contain exactly one molecule definition.

## Hessian files

`read_hessian(...)` expects a whitespace-separated numeric square matrix.

Example:

```text
2.0 3.0
1.0 -2.0
```

The Hessian file itself does not encode energy units. You choose the conversion when you call:

- `calculate(...; wavenumber = wavenumber_kcal)`
- `calculate(...; wavenumber = wavenumber_hartree)`
- `calculate(...; wavenumber = wavenumber_eV)`

or through the CLI `--unit` flag.
