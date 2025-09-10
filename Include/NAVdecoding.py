import numpy as np
from Include.eph_structure_init import eph_structure_init
from Include.ephemeris import ephemeris
from Common.navPartyChk import nav_party_chk
import Common.xcorr as mlxcorr

def NAVdecoding(I_P_InputBits, settings):
    """
    Finds the first preamble occurrence in the bit stream of each channel.
    The preamble is verified by checking the spacing between preambles (6sec)
    and parity checking of the first two words in a subframe. At the same time,
    returns list of channels that are in tracking state and with valid preambles
    in the nav data stream.

    Inputs:
        I_P_InputBits   - output from the tracking function

    Outputs:
        subFrameStart   - Starting position of the first message in the input bit stream.
                          The position is CNAV bit (20ms before convolutional decoding)
                          count since start of tracking. Set to inf if no valid preambles
                          were detected in the channel.
        TOW             - Time Of Week (TOW) of the first message (in seconds).
                          Set to inf if no valid preambles were detected in the channel.
        eph             - SV ephemeris.
    """

    #--------------------------------------------------------------------------
    #                         CU Multi-GNSS SDR  
    # (C) Written by Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
    #--------------------------------------------------------------------------
    # This program is free software; you can redistribute it and/or
    # modify it under the terms of the GNU General Public License
    # as published by the Free Software Foundation; either version 2
    # of the License, or (at your option) any later version.
    #
    # This program is distributed in the hope that it will be useful,
    # but WITHOUT ANY WARRANTY; without even the implied warranty of
    # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    # GNU General Public License for more details.
    #
    # You should have received a copy of the GNU General Public License
    # along with this program; if not, write to the Free Software
    # Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
    # USA.
    #--------------------------------------------------------------------------
    # CVS record:
    # $Id: findPreambles.m,v 1.1.2.10 2017/01/19 21:13:22 dpl Exp $

    #--- Initialize ephemeris structure --------------------------------------
    # This is to ensure variable 'eph' for each SV has a similar structure
    # when only one or even none of the three requisite messages is decoded for a given PRN.
    eph = eph_structure_init()

    # Starting position of the first message in the input bit stream
    subFrameStart = float('inf')

    # TOW of the first message
    TOW = float('inf')

    # Bit and frame synchronization ====================================
    # Preamble search can be delayed to a later point in the tracking results
    # to avoid noise due to tracking loop transients
    searchStartOffset = 0

    #--- Generate the preamble pattern ----------------------------------------
    preamble_bits = np.array([1, -1, -1, -1, 1, -1, 1, 1], dtype=np.int8)

    # "Upsample" the preamble - make 20 values per one bit. The preamble must be
    # found with precision of a sample.
    preamble_ms = np.kron(preamble_bits, np.ones(20, dtype=np.int8))

    # Correlate tracking output with preamble =================================
    # Read output from tracking. It contains the navigation bits. The start
    # of record is skipped here to avoid tracking loop transients.
    bits = np.array(I_P_InputBits[searchStartOffset:], dtype=float)

    # Now threshold the output and convert it to -1 and +1
    bits[bits > 0] = 1
    bits[bits <= 0] = -1

    # Correlate tracking output with the preamble
    tlmXcorrResult, _ = mlxcorr.xcorr(bits, preamble_ms, scale='none')

    # Find all starting points of all preamble-like patterns =================
    xcorrLength = (len(tlmXcorrResult) + 1) // 2

    #--- Find at what index/ms the preambles start ------------------------
    index = np.where(np.abs(tlmXcorrResult[xcorrLength - 1:]) > 153)[0] + searchStartOffset

    # Delete index inferior to 40 and superior to msToProcess-(20*60-1) to avoid boundaries problems
    index = index[(index > 40) & (index < settings.msToProcess - (20 * 60 - 1))]

    # Analyze detected preamble-like patterns ================================
    for i in range(len(index)):
        # Find distances in time between this occurrence and the rest of preamble-like patterns.
        # If the distance is 6000 milliseconds (one subframe), do further verifications by validating the parities of two GPS words
        index2 = index - index[i]
        if np.any(index2 == 6000):
            # Re-read bit values for preamble verification
            # Preamble occurrence is verified by checking the parity of the first two words in the subframe.
            # Now it is assumed that bit boundaries are known. Therefore the bit values over 20ms are combined
            # to increase receiver performance for noisy signals.
            # In total 62 bits must be read:
            # 2 bits from previous subframe are needed for parity checking;
            # 60 bits for the first two 30bit words (TLM and HOW words).
            # The index is pointing at the start of TLM word.
            bits_window = I_P_InputBits[index[i] - 40 : index[i] + 20 * 60]
            bits_window = np.array(bits_window, dtype=float)
            bits_window = bits_window.reshape((20, -1), order='F')
            bits_sum = np.sum(bits_window, axis=0)

            # Now threshold and make it -1 and +1
            bits_sum[bits_sum > 0] = 1
            bits_sum[bits_sum <= 0] = -1

            # Check the parity of the TLM and HOW words
            if nav_party_chk(bits_sum[:32]) != 0 and nav_party_chk(bits_sum[30:62]) != 0:
                # Parity was OK. Record the preamble start position. Skip the rest of preamble pattern checking for this channel
                subFrameStart = index[i]
                break

    # Exclude channel from the active channel list if no valid preamble was detected
    if subFrameStart == float('inf'):
        print('Could not find valid preambles in channel!')
        return eph, subFrameStart, TOW

    # Decode ephemerides ===============================================
    # Convert tracking output to navigation bits =======================
    # Copy 5 sub-frames long record from tracking output ---------------
    navBitsSamples = I_P_InputBits[(subFrameStart - 20) + 1 : (subFrameStart + (1500 * 20)) + 1]
    navBitsSamples = np.array(navBitsSamples, dtype=float)

    # Group every 20 values of bits into columns ------------------------
    navBitsSamples = navBitsSamples.reshape((20, -1), order='F')

    # Sum all samples in the bits to get the best estimate -------------
    navBits = np.sum(navBitsSamples, axis=0)

    # Now threshold and make 1 and 0 -----------------------------------
    # The expression (navBits > 0) returns an array with elements set to 1 if the condition is met and set to 0 if it is not met.
    navBits = (navBits > 0).astype(int)

    # Convert from decimal to binary -----------------------------------
    # The function ephemeris expects input in binary form. In Matlab it is a string array containing only "0" and "1" characters.
    navBitsBin = np.array([str(b) for b in navBits])

    # Decode ephemerides and TOW of the first sub-frame ================
    eph, TOW = ephemeris(navBitsBin[1:1501], navBitsBin[0])
    return eph, subFrameStart, TOW
