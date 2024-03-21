[![codecov](https://codecov.io/gh/MolarVerse/VibrationalAnalysis.jl/graph/badge.svg?token=kESDHEzXcY)](https://codecov.io/gh/MolarVerse/VibrationalAnalysis.jl)
[![CI](https://github.com/MolarVerse/VibrationalAnalysis.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MolarVerse/VibrationalAnalysis.jl/actions/workflows/CI.yml)
[![TagBot](https://github.com/MolarVerse/VibrationalAnalysis.jl/actions/workflows/TagBot.yml/badge.svg)](https://github.com/MolarVerse/VibrationalAnalysis.jl/actions/workflows/TagBot.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10829271.svg)](https://doi.org/10.5281/zenodo.10829271)



# VibrationalAnalysis.jl

This package provides tools to perform [vibrational analysis](https://gaussian.com/vib/) of molecules. It is based on the output of the `QMCFC` molecular dynamics code developed in the [Hofer Lab](https://www.uibk.ac.at/en/aatc/ag-hofer/) of the [University of Innsbruck](https://www.uibk.ac.at/).

## Installation
```julia-repl
julia> ] add VibrationalAnalysis
```

## Usage
```julia-repl
julia> using VibrationalAnalysis
julia> read_calculate("restart.rst", "hessian.dat", "moldescriptor.dat")
julia> calculate(atom_masses, atom_coords, atom_charges, hessian)
```

## Documentation
```julia-repl
julia> using VibrationalAnalysis
julia> ?VibrationalAnalysis
```
