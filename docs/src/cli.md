# Intallation of the CLI

If you want to use the command line interface, you can install it by running the following command:

```julia-repl
julia> using VibrationalAnalysis
julia> VibrationalAnalysis.comonicon_install()
```
Make sure you have `~/.julia/bin` in your PATH.

or by running the following command:

```bash
julia -e 'using Pkg; Pkg.add("VibrationalAnalysis"); using VibrationalAnalysis; VibrationalAnalysis.comonicon_install()'
```

or by build the package from in your default environment:

```julia-repl
julia> using Pkg
julia> Pkg.build("VibrationalAnalysis")
```

## Usage

Run the following command to see the available options:

```bash
vibrationalanalysis -h
```
