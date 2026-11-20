# ------------------------------------------------------------------------------
# Test suite
# ------------------------------------------------------------------------------

using Test, Documenter, Pulvis
DocMeta.setdocmeta!(Pulvis, :DocTestSetup, :(using Pulvis, LinearAlgebra); recursive = true)

# ------------------------------------------------------------------------------
# Documentation tests
# ------------------------------------------------------------------------------

@testset "  documentation tests                        " begin
    doctest(Pulvis; manual = false)
end

# ------------------------------------------------------------------------------
