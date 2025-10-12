import numpy as np
import pandas as pd
from init_settings import init_settings

def generate_L1c_code(PRN):
    """
    generate_L1c_code generates two L1c codes - one for the Pilot signal
    and one for the Data signal.  Each of these is a numpy array that is
    10230 chips long.

    It does _not_ include the TMBOC phase shifts.

    CAcode = generate_ca_code(PRN)

    Inputs:
        PRN         - PRN number of the sequence.
    Outputs:
        L1c_code      - a tuple containing the L1c pilot and data spreading sequences

    --------------------------------------------------------------------------
    This program written by Joe Carey for ASEN 6090 at University of Colorado
    at Boulder.

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
    USA.
    --------------------------------------------------------------------------


    """
    

import numpy as np

def jacobi(n: int, k: int) -> int:
    """
    Computes the Jacobi symbol (n/k) for odd k > 0.
    Returns -1, 0, or +1.
    """
    assert k > 0 and k % 2 == 1
    n = n % k
    t = 1
    while n != 0:
        # factor out powers of 2 from n
        while n % 2 == 0:
            n //= 2
            r = k % 8
            if r == 3 or r == 5:
                t = -t
        # quadratic reciprocity step
        n, k = k, n
        if (n % 4 == 3) and (k % 4 == 3):
            t = -t
        n = n % k
    return t if k == 1 else 0
#
#
#
def jacobi_full_sequence(k: int) -> np.ndarray:
    """
    Compute the Jacobi sequence L(n) for n = 0..k-1.
    Values returned are in {-1, 0, +1}.

    Jacobi is a more general version of the Legendre sequence.
    Legendre sequence results when Jacobi called with a prime number.

    Parameters
    ----------
    k : int
        Odd modulus. For GPS L1C, k = 10223.

    Returns
    -------
    L : np.ndarray
        Array of shape (k,), with entries in {-1,0,+1}.
    """
    values = []
    for n in range(k):
        values.append(jacobi(n, k))  
    return np.array(values, dtype=np.int8)

def weil_sequence(W_index: int, P:int) -> np.ndarray:
    k = 10223
    L_pm1 = jacobi_full_sequence(k).astype(np.int8)   #Legendre sequence returns {-1, 0, +1}
    L_bin = (L_pm1 == 1).astype(np.uint8)             #binary Legendre sequence: map {-1, 0} to 0, +1 to 1
    W1 = W_index % k                                   #handle the edge cases
    L_shift = np.roll (L_bin, -1*W1)
    w = L_bin ^ L_shift
    expansion = np.array([-1, 1, 1, -1, 1, -1, -1])
    break_index = P-1
    if not (0 <= break_index < k) or not (0 < break_index <= k):
        raise ValueError("P out of range for 7-chip replacement.")
    weil = np.concatenate((w[:break_index], expansion, w[break_index:]))
    return (weil)

def bits_to_octal(bits01: np.ndarray, msb_first: bool = True, reverse_bits_in_group: bool = False) -> str:
    """Pack bitstream of 0/1 bits to octal string, 3 bits per digit."""
    # truncate to multiple of 3 for clean comparison
    n = (bits01.size // 3) * 3
    b = bits01[:n].reshape(-1, 3)
    if reverse_bits_in_group:
        b = b[:, ::-1]
    if msb_first:
        vals = (b[:,0] << 2) | (b[:,1] << 1) | b[:,2]
    else:
        # if someone meant LSB-first groups
        vals = b[:,0] | (b[:,1] << 1) | (b[:,2] << 2)
    return "".join(str(int(v)) for v in vals)


def L1c_FS (prn: int, settings) -> np.ndarray:
    """
    Function to produce L1c pilot (including BOC) at the necessary sampling rate 

    Function receives:
        a prn number corresponding to an SV
        a set of settings for the data file & processing
    Returns:
        a numpy array with the L1 pilot signal (including BOC(1,1)
        at the sample rate of the data.  
        The L1c pilot signal will be 10mS long
    """
    df = pd.read_csv("L1Cp.csv",
                 dtype={"PRN": int,
                        "w_pilot": int,
                        "p_pilot": int,
                        "Initial_pilot_code": str,
                        "Final_pilot_code":str
                       }
                )
    L1c_table = df.set_index("PRN", verify_integrity=True)
    W1 = L1c_table.at[PRN, "w_pilot"]
    P1 = L1c_table.at[PRN, "p_pilot"]
    print (f'W1 = {W1}, P1={P1}')
    weil_code = weil_sequence(W1, P1)
    weil_bpsk = 2*ws-1               #re-map the weil code from binary to {+1, -1}
    weil_complement = -1*weil_bpsk   #invert the bpsk version
    weil_combined = np.column_stack((weil_bpsk, weil_complement))
    BOC11 = weil_combined.ravel()
    #
    #  Figure out the stretch factors to re-map the BOC11 to make it 10mS long
    #  At the sampling rate of the captured signal
    #
    codeLength = len(weil_code)
    codePeriod = 0.01   # Period of the L1c code in seconds
    samplingFreq = settings.samplingFreq
    samples_per_code = round((codePeriod * samplingFreq))
    #
    #
    # We use 2*codeLength because now this is a BOC11 - the weil code got twice as long
    index_float = (((2 * codeLength)/ (samples_per_code - 1)) * np.arange(0, samples_per_code) - 0.5)
    index = (np.round(index_float)).astype(int)
    index[-1] = samples_per_code   # just in case round() rolled the last index out of bounds
    L1c_samples = BOC11[index]
    return (L1c_samples)