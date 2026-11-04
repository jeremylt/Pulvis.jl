# ------------------------------------------------------------------------------
# Pulvis types
# ------------------------------------------------------------------------------

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
