# ------------------------------------------------------------------------------
# density
# ------------------------------------------------------------------------------

function computedensity_particle!(
    numneighbors::Int,
    neighbors::Array,
    dimension::Int,
    position::Array,
    radius::Array,
    mass::Array,
    density::Array,
)
    particle = threadIdx().x

    if particle <= numparticles
        sumweightedmassself!(particle, radius, mass, density)
        sumweightedmassneighbors!(
            particle,
            numneighbors[particle],
            neighbors,
            dimension,
            position,
            radius,
            mass,
            density,
        )
    end
end

function computedensity_cuda!(
    numneighbors::Int,
    neighbors::Array,
    dimension::Int,
    position::Array,
    radius::Array,
    mass::Array,
    density::Array,
)
    @cuda threads = length(density) computedensity_particle!(
        numneighbors,
        neighbors,
        dimension,
        position,
        radius,
        mass,
        density,
    )
end

# ------------------------------------------------------------------------------
# pressure
# ------------------------------------------------------------------------------

function computepressure_particle!(
    density,
    initialdensity,
    speedofsound,
    refdensityfactor,
    exponent,
    pressure::Array,
)
    particle = threadIdx().x

    if particle <= numparticles
        densitytopressure!(
            particle,
            density,
            initialdensity,
            speedofsound,
            refdensityfactor,
            exponent,
            pressure,
        )
    end
end

function computepressure_cuda!(
    density,
    initialdensity,
    speedofsound,
    refdensityfactor,
    exponent,
    pressure::Array,
)
    @cuda threads = length(pressure) computepressure_particle!(
        density,
        initialdensity,
        speedofsound,
        refdensityfactor,
        exponent,
        pressure,
    )
end

# ------------------------------------------------------------------------------
# momentum and forces
# ------------------------------------------------------------------------------

function computeforces_cuda!(
    numneighbors::Int,
    neighbors::Array,
    typeids::Array,
    dynamicviscosity::Array,
    bulkviscosity::Array,
    dimension::Int,
    position::Array,
    radius::Array,
    mass::Array,
    density::Array,
    pressure::Array,
    velocity::Array,
    acceleration::Array,
    modifedvelocity::Array,
    modifiedacceleration::Array,
)
    particle = threadIdx().x

    if particle <= numparticles
        momentumequationneighbors!(
            particle,
            numneighbors[particle],
            neighbors,
            typeids,
            dynamicviscosity,
            bulkviscosity,
            position,
            radius,
            mass,
            density,
            pressure,
            velocity,
            acceleration,
            modifedvelocity,
            modifiedacceleration,
        )
    end
end

function computeforces_cuda!(
    numneighbors::Int,
    neighbors::Array,
    typeids::Array,
    dynamicviscosity::Array,
    bulkviscosity::Array,
    position::Array,
    radius::Array,
    mass::Array,
    density::Array,
    pressure::Array,
    velocity::Array,
    acceleration::Array,
    modifedvelocity::Array,
    modifiedacceleration::Array,
)
    @cuda threads = (length(modifiedacceleration)/dimension) accelerationaddmomentumcontribution_particle!(
        numneighbors,
        neighbors,
        typeids,
        dynamicviscosity,
        bulkviscosity,
        position,
        radius,
        mass,
        density,
        pressure,
        velocity,
        acceleration,
        modifedvelocity,
        modifiedacceleration,
    )
end
