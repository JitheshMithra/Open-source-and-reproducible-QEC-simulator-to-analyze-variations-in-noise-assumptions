# decoding for the repetition code, majority vote decoder

def encodeValidation(encodedBit):
    if encodedBit is None or len(encodedBit) == 0:
        raise ValueError("The encoded bit list has to be non empty")
    for value in encodedBit:
        if value not in (0, 1):
            raise ValueError("The encoded bit list must have only 0 and 1")

def majorityDecoder(encodedBit):
    # decode a repetition code with majority vote, but this requires the code length to be odd
    encodeValidation(encodedBit)
    n = len(encodedBit)
    # this will require odd n to avoid ties
    if n % 2 == 0:
        raise ValueError("Repetition code length must be odd")
    oneCounter = 0
    for value in encodedBit:
        if value == 1:
            oneCounter = oneCounter + 1
    if oneCounter > n // 2:
        return 1
    else:
        return 0
