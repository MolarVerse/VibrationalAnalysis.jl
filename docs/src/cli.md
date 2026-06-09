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

Run a calculation with restart and Hessian files:

```bash
vibrationalanalysis restart.rst hessian.dat --unit kcal --output wavenumbers.dat
```

Run the same calculation from a standard single-structure XYZ file:

```bash
vibrationalanalysis structure.xyz hessian.dat --unit kcal --output wavenumbers.dat
```

Include intensities by passing a moldescriptor file:

```bash
vibrationalanalysis restart.rst hessian.dat --moldescriptor moldescriptor.dat --output wavenumbers.dat
```

When an XYZ file is used, atom types default to `1` for all atoms. Multi-molecule moldescriptor files therefore still require restart input with explicit atom-type information.

Write normal modes in matrix notation or xyz trajectories:

```bash
vibrationalanalysis restart.rst hessian.dat --normal-modes normal_modes.dat --modes
```
