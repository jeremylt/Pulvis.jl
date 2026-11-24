# ------------------------------------------------------------------------------
# Pulvis - Automatic Differentiated, GPU Enabled Particle Physics
# ------------------------------------------------------------------------------

module Pulvis

# ------------------------------------------------------------------------------
# Standard libraries
# ------------------------------------------------------------------------------

using Base.Threads
using CUDA
using Adapt
using LinearAlgebra: norm, dot, ⋅
using StaticArrays

# ------------------------------------------------------------------------------
# Imports
# ------------------------------------------------------------------------------

# Types
include("Types.jl")

# Utilities
include("Binning.jl")
include("Data.jl")
include("Enums.jl")
include("Kernels.jl")
include("OpenBoundary.jl")
include("Setup.jl")

# Physics
include("EquationsOfState.jl")
include("Density.jl")
include("MomentumFormulations.jl")
include("Momentum.jl")
include("InteractionModels.jl")

# Time stepping
include("TimeStep.jl")

# Drivers
include("Drivers.jl")

# ------------------------------------------------------------------------------
# User available types and methods
# ------------------------------------------------------------------------------

export PulvisScalar, StateData, PulvisVector
export TypeInfo, StateInfo, BinInfo
export Kernel, MomentumFormulation, EquationOfState, InteractionModel
export setkernelshape, setmomentumformulation, setequationsofstate
export setdevicecomputation

end # module Pulvis
