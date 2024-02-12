""" 
# VibrationalAnalysis.jl 

This module contains functions to perform vibrational analysis on a QMCFC output.

## Usage

Can read directly from rst-file, moldescriptor-file and hessian-file and perform vibrational analysis on the system.

```julia-repl
julia> read_calculate("restart.rst", "hessian.dat", "moldescriptor.dat")
```

Or directly perform a vibrational analysis with atom masses, atom coordinates, atom charges and hessian of the system.

```julia-repl
julia> read_calculate(atom_masses, atom_coords, atom_charges, hessian)
```

"""
module VibrationalAnalysis

using LinearAlgebra

include("read_files.jl")
include("coordinate_transform.jl")
include("symmetrize.jl")
include("transformation.jl")
include("observables.jl")
include("write_out.jl")
include("calculate.jl")

end # module VibrationalAnalysis
