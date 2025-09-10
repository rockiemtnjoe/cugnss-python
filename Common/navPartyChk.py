import numpy as np

def nav_party_chk(ndat):
    """
    This function computes and checks the parity bits on a GPS navigation word.
    Based on the flowchart in Figure 2-10 in the 2nd Edition of the GPS-SPS Signal Spec.

    Parameters
    ----------
    ndat : array_like
        Array (length 32) of 32 bits representing a GPS navigation word:
        bits [-2, -1, 0, 1, ..., 28, 29] (D29*, D30*, d1-d24, D25-D30).

    Returns
    -------
    int
        Test value: +1 or -1 if parity PASSED, 0 if parity fails.
        +1 means bits #1-24 have correct polarity, -1 means bits #1-24 must be inverted.
    """

    # In order to accomplish the exclusive or operation using multiplication,
    # this program represents a '0' with a '-1' and a '1' with a '1'.
    # See the MATLAB comments for truth table.

    ndat = np.array(ndat, dtype=int)

    #--- Check if the data bits must be inverted ------------------------------
    # If D30* (ndat[1]) is not 1, invert bits d1-d24 (ndat[2:26])
    if ndat[1] != 1:
        ndat[2:26] = -1 * ndat[2:26]

    #--- Calculate 6 parity bits ----------------------------------------------
    # The elements of ndat correspond to the bits in ICD-200C Table 20-XIV:
    # ndat[0] = D29*, ndat[1] = D30*, ndat[2:26] = d1-d24, ndat[26:32] = received D25-D30
    # parity contains computed D25-D30 bits

    parity = np.zeros(6, dtype=int)
    parity[0] = ndat[0] * ndat[2] * ndat[3] * ndat[4] * ndat[6] * ndat[7] * ndat[11] * ndat[12] * ndat[13] * ndat[14] * ndat[15] * ndat[18] * ndat[19] * ndat[21] * ndat[24]
    parity[1] = ndat[1] * ndat[3] * ndat[4] * ndat[5] * ndat[7] * ndat[8] * ndat[12] * ndat[13] * ndat[14] * ndat[15] * ndat[16] * ndat[19] * ndat[20] * ndat[22] * ndat[25]
    parity[2] = ndat[0] * ndat[2] * ndat[4] * ndat[5] * ndat[6] * ndat[8] * ndat[9] * ndat[13] * ndat[14] * ndat[15] * ndat[16] * ndat[17] * ndat[20] * ndat[21] * ndat[23]
    parity[3] = ndat[1] * ndat[3] * ndat[5] * ndat[6] * ndat[7] * ndat[9] * ndat[10] * ndat[14] * ndat[15] * ndat[16] * ndat[17] * ndat[18] * ndat[21] * ndat[22] * ndat[24]
    parity[4] = ndat[1] * ndat[2] * ndat[4] * ndat[6] * ndat[7] * ndat[8] * ndat[10] * ndat[11] * ndat[15] * ndat[16] * ndat[17] * ndat[18] * ndat[19] * ndat[22] * ndat[23] * ndat[25]
    parity[5] = ndat[0] * ndat[4] * ndat[6] * ndat[7] * ndat[9] * ndat[10] * ndat[11] * ndat[12] * ndat[14] * ndat[16] * ndat[20] * ndat[23] * ndat[24] * ndat[25]

    #--- Compare if the received parity is equal to the calculated parity ------
    # If all parity bits match, return -1 * ndat[1] (D30*), else return 0
    if np.sum(parity == ndat[26:32]) == 6:
        return -1 * ndat[1]
    else:
        return 0
