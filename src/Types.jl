# ------------------------------------------------------------------------------
# Pulvis types
# ------------------------------------------------------------------------------

PulvisScalar = Float64
StateData = SVector
PulvisVector = CuVector

struct TypeInfo
    name::String
    numparticles::Int
end

struct StateInfo
    name::String
    dimension::Int
end

struct BinInfo
    lower::Array{PulvisScalar}
    upper::Array{PulvisScalar}
    width::PulvisScalar
    counts::Array{Int}
end

struct TimeInfo
    step::Int
    current::PulvisScalar
    dt::PulvisScalar
    final::PulvisScalar
end
