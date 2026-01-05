#simple repetition code simulation, for now at least
#estimate logical error rate vs p for majority vote repetition code
import random
from .noise import bitflipnoise, validationp
from .decode import majoritydecoder
def encoderepetition(logicalbit, n):
    #logicalbit has to be 0 or 1, n has to be positive and odd
    if logicalbit not in (0,1):
        raise ValueError("Logical bit must be 0 or 1")
    if n <= 0:
        raise ValueError("n must be positive")
    if n%2==0:
        raise ValueError("n must be odd")
    return [logicalbit]*n

def runtrials(n, p, trials, seed, logicalbit):
    #validating the inputs
    if trials <= 0:
        raise ValueError("trials must be positive")
    validationp(p)
    rng=random.Random(seed)
    failures=0
    for i in range(trials):
        encodedbit=encoderepetition(logicalbit, n)
        noisedbit= bitflipnoise(encodedbit, p, rng)
        decodedbit= majoritydecoder(noisedbit)
        if decodedbit!=logicalbit:
            failures=failures+1
    return failures/trials

def psweep(n, pvalues, trials, seed, logicalbit):
    results=[]
    for i in pvalues:
        logicalerrorrate= runtrials(n, i, trials, seed + int(i*1e6), logicalbit)
        results.append((i, logicalerrorrate))
    return results