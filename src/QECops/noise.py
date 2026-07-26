# noise models for the repetition code simulator

def validation(encodedBit):
    # encodedBit is and has to be a non empty list, a sequence of 0 and 1
    if encodedBit is None or len(encodedBit) == 0:
        raise ValueError("The encoded bit list has to be non empty")
    for value in encodedBit:
        if value not in (0, 1):
            raise ValueError("The encoded bit list must have only 0 and 1")

def validationP(p):
    if p < 0 or p > 1:
        raise ValueError("Probability/p must be between 0 and 1")

def bitflipNoise(bits, p, rng):
    validation(bits)
    validationP(p)
    out = []
    for b in bits:
        if rng.random() < p:
            out.append(1 - b)
        else:
            out.append(b)
    return out

def depolarizingNoise(bits, p, rng):
    # right now symmetric flip, but will be extended to x/z channels if needed when quantum states come along in later versions
    validation(bits)
    validationP(p)

    effectiveP = 2 * p / 3
    return bitflipNoise(bits, effectiveP, rng)

def biasedNoise(bits, px, pz, rng):
    validation(bits)
    validationP(px)
    validationP(pz)
    if px + pz > 1:
        raise ValueError("px + pz must be <= 1")

    out = []
    for b in bits:
        r = rng.random()
        if r < px:
            out.append(1 - b)  # bit flip
        # Z errors do not affect classical bit values in repetition code
        elif r < px + pz:
            out.append(b)  # phase flip (no change in bit value)
        else:
            out.append(b)  # no error

    return out

def correlatedNoise(bits, p, correlation, rng):
    validation(bits)
    validationP(p)
    validationP(correlation)
    n = len(bits)
    newBits = bits.copy()

    i = 0
    while i < n:
        if rng.random() < p:
            newBits[i] = 1 - newBits[i]  # flip the bit

            # propagate the correlation
            j = i + 1
            while j < n and rng.random() < correlation:
                newBits[j] = 1 - newBits[j]
                j += 1
            i = j
        else:
            i += 1
    return newBits

def applyNoise(bits, rng, noiseType="depolarizing", **params):
    # dispatches to the right noise function based on noiseType, params depend on which model is picked (p for bitflip/depolarizing, px/pz for biased, p+correlation for correlated)
    if noiseType == "bitflip":
        return bitflipNoise(bits, params["p"], rng)
    elif noiseType == "depolarizing":
        return depolarizingNoise(bits, params["p"], rng)
    elif noiseType == "biased":
        return biasedNoise(bits, params["px"], params["pz"], rng)
    elif noiseType == "correlated":
        return correlatedNoise(bits, params["p"], params["correlation"], rng)
    else:
        raise ValueError(f"Unknown noise type: {noiseType}")
