# ------------------------------------------------------------------------------
# Pulvis types
# ------------------------------------------------------------------------------

PulvisScalar = Float64
StateData = SVector
PulvisVector = CuVector

struct TypeInfo
    name::String
    count::Int
end

struct StateInfo
    name::String
    dimension::Int
end

struct BinInfo
    lower::Array{PulvisScalar}
    upper::Array{PulvisScalar}
    binwidth::PulvisScalar
    numbins::Array{Int}
end
