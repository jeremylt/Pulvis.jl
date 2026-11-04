# ------------------------------------------------------------------------------
# Pulvis - Automatic Differentiated, GPU Enabled Particle Physics
# ------------------------------------------------------------------------------

module Pulvis

# Dependencies
using CUDA
using Adapt
using LinearAlgebra: norm, dot, ⋅
using StaticArrays

# Types
include("Types.jl")

# Utilities
include("Data.jl")
include("Enums.jl")
include("Kernels.jl")

# Physics
include("EquationsOfState.jl")
include("Density.jl")
include("MomentumFormulations.jl")
include("Momentum.jl")

# Drivers
include("Drivers.jl")

end # module Pulvis
