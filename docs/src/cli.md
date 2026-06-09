# Installation of the CLI

If you want to use the command line interface, you can install it by running the following command:

```julia-repl
julia> using VibrationalAnalysis
julia> VibrationalAnalysis.comonicon_install()
```
Make sure you have `~/.julia/bin` in your PATH.

You can also install the package and CLI in one shell command:

```bash
julia -e 'using Pkg; Pkg.add("VibrationalAnalysis"); using VibrationalAnalysis; VibrationalAnalysis.comonicon_install()'
```

Or build the package from your default environment:

```julia-repl
julia> using Pkg
julia> Pkg.build("VibrationalAnalysis")
```

## Usage

Run the following command to see the available options:

```bash
vibrationalanalysis -h
```

## Common runs

Run a calculation with structure and Hessian files:

```bash
vibrationalanalysis structure.rst hessian.dat --unit kcal --output wavenumbers.dat
```

Run the same calculation from a standard single-structure XYZ file:

```bash
vibrationalanalysis structure.xyz hessian.dat --unit kcal --output wavenumbers.dat
```

Include intensities by passing a moldescriptor file:

```bash
vibrationalanalysis structure.rst hessian.dat --moldescriptor moldescriptor.dat --output wavenumbers.dat
```

When an XYZ file is used, atom types default to `1` for all atoms. If `--moldescriptor` is used with an XYZ file, the moldescriptor file must contain exactly one molecule definition because XYZ files do not encode molecule-type assignments.

Write the normal-mode matrix to a file:

```bash
vibrationalanalysis structure.rst hessian.dat --normal-modes normal_modes.dat
```

Write normal modes in matrix notation or xyz trajectories:

```bash
vibrationalanalysis structure.rst hessian.dat --normal-modes normal_modes.dat --modes
```

The `--modes` flag writes XYZ trajectories named `modes-1.xyz`, `modes-2.xyz`, and so on in the current working directory.

## Output files

Without `--moldescriptor`, the main output has three columns:

```text
# Wavenumbers (cm-1)  Force constants (mdyn Å-1)  Reduced masses (amu)
```

With `--moldescriptor`, intensities are added:

```text
# Wavenumbers (cm-1)  Intensities (km mol-1)  Force constants (mdyn Å-1)  Reduced masses (amu)
```

If `--output` is omitted, the table is written to stdout.

## Units

The CLI supports three Hessian unit conventions:

- `--unit kcal`
- `--unit hartree`
- `--unit eV`

The default is `kcal`.

## Common failure cases

These are the most important current error cases to know upfront:

- `The xyz file contains more than one structure.`
  Use a single-structure XYZ file, not a trajectory or multi-frame file.

- `XYZ input requires a moldescriptor file with exactly one molecule type because XYZ files do not encode molecule-type assignments.`
  Use a restart file when molecule-type assignments matter, or reduce the moldescriptor file to one molecule block.

- `Invalid unit. Options are kcal, hartree and eV.`
  The CLI only accepts those three spellings.
