# density

function computedensity!(numneighbors, neighbors, position, radius, mass, density)
    particle = threadIdx().x

    if particle <= numparticles
        sumweightedmassself!(particle, radius, mass, density)
        sumweightedmassneighbors!(particle, numneighbors[particle], neighbors, position, radius, mass, density)
    end
end

# pressure

function computepressure!()
    particle = threadIdx().x

    if particle <= numparticles
        densitytopressure!(particle, density, initialdensity, view(pressure, particle:particle+1))
    end
end

# momentum

function accelerationaddmomentumcontribution!(numneighbors, neighbors, typeids, dynamicviscosity, bulkviscosity, position, radius, mass, density, pressure, velocity, acceleration, modifedvelocity, modifiedacceleration)
    particle = threadIdx().x

    if particle <= numparticles
        momentumequationneighbors!(particle, numneighbors[particle], neighbors, typeids, dynamicviscosity, bulkviscosity, position, radius, mass, density, pressure, velocity, acceleration, modifedvelocity, modifiedacceleration)
    end    
end
