# Migration to 0.4

Version `0.4.0` is a structure-input release. The main changes are about naming and input semantics, not a rewrite of the core vibrational analysis routines.

## What changed

- `read_structure(...)` is now the recommended public structure reader
- the CLI is documented in terms of `structure` input, not restart-only input
- XYZ handling is explicit and documented
- XYZ plus multi-molecule moldescriptor files now raises a clear error instead of relying on ambiguous atom typing

## Recommended code updates

Prefer this:

```julia
atom_names, atom_masses, atom_coords, atom_types = read_structure("structure.rst")
atom_names, atom_masses, atom_coords, atom_types = read_structure("structure.xyz")
```

instead of writing new code directly against format-specific readers unless you genuinely need that control.

## What did not change

- `read_rst(...)` still works
- `read_xyz(...)` still works
- existing restart-based workflows remain valid
- the `calculate(...)` and `write_*` APIs remain the same

## When to keep using restart files

Prefer restart files when:

- your workflow depends on molecule-type assignments
- your moldescriptor file contains multiple molecule definitions
- you want the structure source to carry `atom_types`

## When XYZ is enough

XYZ input is a good fit when:

- you only need geometry and masses
- intensities are not needed
- intensities are needed but the moldescriptor file contains exactly one molecule definition
