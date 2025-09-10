import numpy as np
from .generateCAcode import generate_ca_code

def make_ca_table(PRN, settings):
    """
    Generates digitized C/A codes for the specified PRN based on receiver settings.
    The codes are digitized at the sampling frequency specified in the settings object.
    One array in the returned ca_codes_table is one C/A code. The array index corresponds
    to the PRN number of the C/A code.

    Inputs:
        PRN             - specified PRN for C/A code
        settings        - receiver settings object with attributes:
                            samplingFreq     - sampling frequency (Hz)
                            codeFreqBasis    - code frequency basis (Hz)
                            codeLength       - length of C/A code (chips)
    Outputs:
        ca_codes_table  - numpy array containing digitized C/A code for the given PRN

    --------------------------------------------------------------------------
    SoftGNSS v3.0

    Copyright (C) Darius Plausinaitis
    Written by Darius Plausinaitis
    Based on Peter Rinder and Nicolaj Bertelsen
    --------------------------------------------------------------------------
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

    CVS record:
    $Id: makeCaTable.m,v 1.1.2.6 2006/08/14 11:38:22 dpl Exp $
    """

    #--- Find number of samples per spreading code ----------------------------
    samples_per_code = int(round(
        settings.samplingFreq / (settings.codeFreqBasis / settings.codeLength)
    ))

    #--- Find time constants --------------------------------------------------
    ts = 1 / settings.samplingFreq   # Sampling period in sec
    tc = 1 / settings.codeFreqBasis  # C/A chip period in sec

    #--- Generate CA code for given PRN ---------------------------------------
    ca_code = generate_ca_code(PRN)

    #=== Digitizing ===========================================================
    #--- Make index array to read C/A code values -----------------------------
    # The length of the index array depends on the sampling frequency -
    # number of samples per millisecond (because one C/A code period is one
    # millisecond).
    code_value_index = np.ceil(
        ts * np.arange(1, samples_per_code + 1) / tc
    ).astype(int)

    #--- Correct the last index (due to number rounding issues) ---------------
    code_value_index[-1] = 1023

    #--- Make the digitized version of the C/A code ---------------------------
    # The "upsampled" code is made by selecting values from the CA code
    # chip array (ca_code) for the time instances of each sample.
    ca_codes_table = ca_code[code_value_index - 1]

    return ca_codes_table
