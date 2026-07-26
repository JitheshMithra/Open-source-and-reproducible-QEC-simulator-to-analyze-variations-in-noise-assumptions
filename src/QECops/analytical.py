# engine for producing analytical binomial distribution.
from math import comb

def analyticalLogicalError(n, p):
    if n <= 0:
        raise ValueError("n must be positive")
    if n % 2 == 0:
        raise ValueError("n must be odd")
    if p < 0 or p > 1:
        raise ValueError("p must be between 0 and 1")

    failuresNeeded = (n // 2) + 1
    total = 0

    for k in range(failuresNeeded, n + 1):
        ways = comb(n, k)
        prob = (p**k) * ((1 - p)**(n - k))
        total = total + ways * prob

    return total
