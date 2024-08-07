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
# Test Calculate Functions
include("calculate.jl")
# Test Write Functions
include("write.jl")
# Test Check Functions
include("check.jl")