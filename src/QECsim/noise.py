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

def bitflipnoise(encodedbit, p, rng):
    validation(encodedbit)
    validationp(p)
    noisedbit=[]
    for value in encodedbit:
        if rng.random() < p:
            noisedbit.append(1-value)
        else:
            noisedbit.append(value)
    return noisedbit