<center><img src="assets/logo.png" width="250"></center>

[![codecov](https://codecov.io/gh/MolarVerse/VibrationalAnalysis.jl/graph/badge.svg?token=kESDHEzXcY)](https://codecov.io/gh/MolarVerse/VibrationalAnalysis.jl)
[![CI](https://github.com/MolarVerse/VibrationalAnalysis.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/MolarVerse/VibrationalAnalysis.jl/actions/workflows/CI.yml)
[![TagBot](https://github.com/MolarVerse/VibrationalAnalysis.jl/actions/workflows/TagBot.yml/badge.svg)](https://github.com/MolarVerse/VibrationalAnalysis.jl/actions/workflows/TagBot.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10829271.svg)](https://doi.org/10.5281/zenodo.10829271)


# VibrationalAnalysis.jl

This package provides tools to perform [vibrational analysis](https://gaussian.com/vib/) of molecules. It is based on the output of QMCFC and [PQ](https://github.com/MolarVerse/PQ) molecular dynamics codes developed of the [University of Innsbruck](https://www.uibk.ac.at/).

## Installation
```
using Pkg
Pkg.add("VibrationalAnalysis")
```

## Usage
```julia-repl
julia> using VibrationalAnalysis
julia> read_calculate("restart.rst", "hessian.dat", "moldescriptor.dat")
julia> calculate(atom_masses, atom_coords, atom_charges, hessian)
```

# Acknowledgements
This package was developed as part of the [MolarVerse](https://github.com/MolarVerse) project. Significant contributions were made by:
- Josef M. Gallmetzer @galjos
- Jakob Gamper @97gamjak

## Citation
If you use this package in your research, please cite it using the following DOI: [10.5281/zenodo.10829271](https://doi.org/10.5281/zenodo.10829271)
