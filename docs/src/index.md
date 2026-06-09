```@raw html
<img src="assets/logo.png" width="250">
```

[![codecov](https://codecov.io/gh/MolarVerse/VibrationalAnalysis.jl/graph/badge.svg?token=kESDHEzXcY)](https://codecov.io/gh/MolarVerse/VibrationalAnalysis.jl)
[![CI](https://github.com/MolarVerse/VibrationalAnalysis.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MolarVerse/VibrationalAnalysis.jl/actions/workflows/CI.yml)
[![TagBot](https://github.com/MolarVerse/VibrationalAnalysis.jl/actions/workflows/TagBot.yml/badge.svg)](https://github.com/MolarVerse/VibrationalAnalysis.jl/actions/workflows/TagBot.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10829271.svg)](https://doi.org/10.5281/zenodo.10829271)


# VibrationalAnalysis.jl

This package provides tools to perform [vibrational analysis](https://gaussian.com/vib/) of molecules. It is based on the output of QMCFC and [PQ](https://github.com/MolarVerse/PQ) molecular dynamics codes developed at the [University of Innsbruck](https://www.uibk.ac.at/).

Version `0.4.0` accepts either restart files or standard single-structure XYZ files through the public `read_structure(...)` entrypoint and the CLI.

## Installation
To install the package, run the following command in the Julia REPL:
```julia-repl
julia> ]
pkg> add VibrationalAnalysis
```

## Installation of the CLI
If you want to use the command line interface, you can install it by running the following command:
```julia-repl
julia> using VibrationalAnalysis
julia> VibrationalAnalysis.comonicon_install()
```
Make sure you have `~/.julia/bin` in your PATH.

## Usage
```julia-repl
shell> vibrationalanalysis -h
```

## Documentation map

- [Quickstart](quickstart.md): end-to-end restart and XYZ workflows with real fixture outputs
- [Input Files](input-files.md): exact file expectations for restart, XYZ, Hessian, and moldescriptor inputs
- [CLI](cli.md): command-line recipes, output files, and common failure cases
- [Examples](examples.md): API-level examples for units, Hessian sign conventions, and writing modes
- [Validation](validation.md): current regression coverage, reference values, and linear/near-linear notes
- [Migration](migration.md): what changed in `0.4.0`

# Acknowledgements
This package was developed as part of the [MolarVerse](https://github.com/MolarVerse) organization. Significant contributions were made by:
- Josef M. Gallmetzer @galjos
- Jakob Gamper @97gamjak
- Thomas S. Hofer

## Citation
If you use this package in your research, please cite it using the following DOI: [10.5281/zenodo.10829271](https://doi.org/10.5281/zenodo.10829271)
