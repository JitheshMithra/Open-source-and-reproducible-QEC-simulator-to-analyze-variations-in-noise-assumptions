#decoding for the repetition code, majority vote decoder

def encodevalidation(encodedbit):
    if encodedbit is None or len(encodedbit) == 0:
        raise ValueError("The encoded bit list has to be non empty")
    for value in encodedbit:
        if value not in (0, 1):
            raise ValueError("The encoded bit list must have only 0 and 1")

def majoritydecoder(encodedbit):
    #decode a repetition code with majority vote, but this requires the code length to be odd
    encodevalidation(encodedbit)
    n= len(encodedbit)
    #this will require odd n to avoid ties
    if n%2==0:
        raise ValueError("Repetition code length must be odd")
    onecounter=0
    for value in encodedbit:
        if value==1:
            onecounter= onecounter+1
    if onecounter > n//2:
        return 1
    else:
        return 0
