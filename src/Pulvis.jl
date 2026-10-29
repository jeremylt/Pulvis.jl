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
include("EquationsOfState.jl")
include("Density.jl")
include("MomentumFormulations.jl")
include("Momentum.jl")

# Drivers
include("Drivers.jl")

end # module Pulvis
