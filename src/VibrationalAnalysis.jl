"""@meta
CurrentModule = VibrationalAnalysis
"""

"""
# VibrationalAnalysis.jl 

This module contains functions to perform vibrational analysis on a QMCFC output.

## Usage

Can read directly from rst-file, moldescriptor-file and hessian-file and perform vibrational analysis on the system.

```julia-repl
julia> atom_names, atom_masses, atom_coords, atom_types = read_rst("restart.rst")
julia> hessian = read_hessian("hessian.dat")
julia> atom_charges = read_moldescriptor("moldescriptor.dat", atom_names, atom_types)
julia> calculate(atom_masses, atom_coords, atom_charges, hessian)
```

Or directly perform a vibrational analysis with atom masses, atom coordinates, atom charges and hessian of the system.

```julia-repl
julia> wavenumbers, intensities, _, _ = calculate(atom_masses, atom_coords, atom_charges, hessian)
julia> write_wavenumber_intensity(wavenumbers, intensities, filename="wavenumbers.dat")
```

"""
module VibrationalAnalysis

using Printf
using LinearAlgebra
using Comonicon

include("massesdict.jl")
include("read_files.jl")
include("coordinate_transform.jl")
include("symmetrize.jl")
include("transformation.jl")
include("observables.jl")
include("write.jl")
include("calculate.jl")
include("check.jl")
include("main.jl")

end # module VibrationalAnalysis
