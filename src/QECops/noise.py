#noise models for the repetition code simulator

def validation(encodedbit):
    #encodedbit is and has to be a non empty list, a sequence of 0 and 1
    if encodedbit is None or len(encodedbit) == 0:
        raise ValueError("The encoded bit list has to be non empty")
    for value in encodedbit:
        if value not in (0, 1):
            raise ValueError("The encoded bit list must have only 0 and 1")

def validationp(p):
    if p<0 or p>1:
        raise ValueError("Probability/p must be between 0 and 1")

def bitflipnoise(bits, p,rng):
    validation(bits)
    validationp(p)
    out = []
    for b in bits:
        if rng.random() < p:
            out.append(1 - b)
        else:
            out.append(b)
    return out

def depolarizingnoise(bits,p,rng):
    #right now symmetric flip, but will be extended to x/z channels if needed when quantum states come along in later versions
    validation(bits)
    validationp(p)
   
    effective_p = 2 * p / 3
    return bitflipnoise(bits, effective_p, rng)

def biasednoise(bits, px, pz, rng):
    validation(bits)
    validationp(px)
    validationp(pz)
    if px + pz > 1:
        raise ValueError("px + pz must be <= 1")
   
    out = []
    for b in bits:
        r = rng.random()
        if r < px:
            out.append(1 - b)  #bit flip
        #Z errors do not affect classical bit values in repetition code
        elif r < px + pz:
            out.append(b)  #phase flip (no change in bit value)
        else:
            out.append(b)  #no error

    return out

def correlatednoise(bits,p,correlation, rng):
    validation(bits)
    validationp(p)
    validationp(correlation)
    n = len(bits)
    newbits=bits.copy()
   
    i=0
    while i<n:
        if rng.random()<p:
            newbits[i]=1-newbits[i] #flip the bit
           
            #propogate the correlation
            j=i+1
            while j<n and rng.random()<correlation:
                newbits[j]=1-newbits[j]
                j +=1
            i=j
        else:
            i +=1
    return newbits

def applynoise(bits, rng, noisetype="depolarizing", **params):
    if noisetype == "bitflip":
        return bitflipnoise(bits, params["p"], rng)
    elif noisetype == "depolarizing":
        return depolarizingnoise(bits, params["p"], rng)
    elif noisetype == "biased":
        return biasednoise(bits, params["px"], params["pz"], rng)
    elif noisetype == "correlated":
        return correlatednoise(bits, params["p"], params["correlation"], rng)
    else:
        raise ValueError(f"Unknown noise type: {noisetype}")