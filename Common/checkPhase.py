def checkPhase(word, D30Star):
    """
    Checks the parity of the supplied 30-bit word and inverts the first 24 bits if needed.
    Args:
        word (list or str): 30-bit word as a list or string of '0'/'1'
        D30Star (str): The last bit of the previous word ('0' or '1')
    Returns:
        list: Word with corrected polarity of the data bits (as a list of '0'/'1')
    """
    if isinstance(word, str):
        word = list(word)
    if D30Star == '1':
        # Invert first 24 bits
        word[:24] = ['0' if b == '1' else '1' for b in word[:24]]
    return word
