[![codecov](https://codecov.io/gh/galjos/VibrationalAnalysis.jl/graph/badge.svg?token=PIG1D1QIEE)](https://codecov.io/gh/galjos/VibrationalAnalysis.jl)
[![Build CI](https://github.com/galjos/VibrationalAnalysis.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/galjos/VibrationalAnalysis.jl/actions/workflows/CI.yml)

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
