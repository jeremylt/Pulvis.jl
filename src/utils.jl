function distancebetweenparticles(positioni, positionj)
    if dimension == 1
        return abs(positionj[1] - positioni[1])
    elseif dimension == 2
        return sqrt((positionj[1] - positioni[1])^2 + (positionj[2] - positioni[2])^2)
    elseif dimension == 3
        return sqrt(
            (positionj[1] - positioni[1])^2 +
            (positionj[2] - positioni[2])^2 +
            (positionj[3] - positioni[3])^3,
        )
    end
end
