# ------------------------------------------------------------------------------
# Pulvis - Automatic Differentiated, GPU Enabled Particle Physics
# ------------------------------------------------------------------------------

module Pulvis

# Dependencies
using LinearAlgebra: norm, dot

# Utilities
include("Enums.jl")
include("Kernels.jl")
include("Utilities.jl")

# Physics
include("Density.jl")
include("EquationsOfState.jl")
include("MomentumFormulations.jl")
include("Momentum.jl")

# Hardware specific implementations
include("CPU.jl")
include("CUDA.jl") # ToDo allow optional usage of CUDA

end # module Pulvis
