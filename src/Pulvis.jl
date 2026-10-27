# ------------------------------------------------------------------------------
# Pulvis - Automatic Differentiated, GPU Enabled Particle Physics
# ------------------------------------------------------------------------------

module Pulvis

# Utilities
include("Enums.jl")
include("Kernel.jl")
include("EquationsOfState.jl")
include("Utilities.jl")

# Physics
include("Density.jl")
include("Momentum.jl")

# Hardware specific implementations
include("CPU.jl")
include("CUDA.jl") # ToDo allow optional usage of CUDA

end # module Pulvis
