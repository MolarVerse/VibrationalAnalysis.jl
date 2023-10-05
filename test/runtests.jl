using VibrationalAnalysis
using Test

# Test reading RST files src/read_files.jl
include("read_files.jl")
# Test Symmetric Matrix Functions
include("symmetrize.jl")
# Test Coordinate Transformation Functions
include("coordinate_transform.jl")
# Test Transformation Functions
include("transformation.jl")
# Test Observable Functions
include("observables.jl")