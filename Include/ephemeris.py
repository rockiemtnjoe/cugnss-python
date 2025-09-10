def ephemeris(bits, D30Star):
    """
    Decodes ephemerides and TOW from the given bit stream (5 subframes, 1500 bits).
    The first element in bits must be the first bit of a subframe.
    The subframe ID of the first subframe is not important.

    Function does not check parity!

    Args:
        bits (list[str] or np.ndarray): Bits of the navigation messages (5 subframes).
            Must contain only '0' or '1' strings, length 1500.
        D30Star (str): The last bit of the previous nav-word ('0' or '1').
            Refer to GPS ICD (IS-GPS-200D) for parity checking details.

    Returns:
        tuple: (eph, TOW)
            TOW (int): Time Of Week (TOW) of the first sub-frame in the bit stream (seconds)
            eph (dict): SV ephemeris

    --------------------------------------------------------------------------
    SoftGNSS v3.0

    Copyright (C) Darius Plausinaitis and Kristin Larson
    Written by Darius Plausinaitis and Kristin Larson
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
    """

    # Check if there is enough data ==========================================
    import numpy as np
    from Include.eph_structure_init import eph_structure_init
    from Common.twosComp2dec import twosComp2dec
    from Common.checkPhase import checkPhase

    if len(bits) < 1500:
        raise ValueError('The parameter bits must contain 1500 bits!')
    if not (isinstance(bits, (list, np.ndarray)) and all(b in {'0', '1'} for b in bits)):
        raise ValueError('The parameter bits must be a list/array of \'0\' and \'1\' strings!')
    if not (isinstance(D30Star, str) and D30Star in {'0', '1'}):
        raise ValueError('The parameter D30Star must be a string "0" or "1"!')

    # Pi used in the GPS coordinate system
    gpsPi = 3.1415926535898

    #--- Initialize ephemeris structure  --------------------------------------
    # Ensures variable 'eph' for each SV has a similar structure even if only one or none of the three requisite messages is decoded.
    eph = eph_structure_init()
    subframe = None
    D30Star_local = D30Star

    # Decode all 5 sub-frames ================================================
    for i in range(5):
        #--- "Cut" one sub-frame's bits ---------------------------------------
        subframe = bits[300*i : 300*(i+1)]
        subframe = list(subframe)

        #--- Correct polarity of the data bits in all 10 words ----------------
        for j in range(10):
            idx0 = 30*j
            idx1 = 30*(j+1)
            subframe[idx0:idx1] = checkPhase(subframe[idx0:idx1], D30Star_local)
            D30Star_local = subframe[idx1-1]

        #--- Decode the sub-frame id ------------------------------------------
        # For more details on sub-frame contents please refer to GPS IS.
        subframeID = int(''.join(subframe[49:52]), 2)

        #--- Decode sub-frame based on the sub-frames id ----------------------
        # Select the necessary bits and convert them to decimal numbers.
        # For more details on sub-frame contents please refer to GPS ICD (IS-GPS-200D).
        if subframeID == 1:
            # It contains WN, SV clock corrections, health and accuracy
            eph['weekNumber'] = int(''.join(subframe[60:70]), 2) + 1024
            eph['accuracy'] = int(''.join(subframe[72:76]), 2)
            eph['health'] = int(''.join(subframe[76:82]), 2)
            eph['T_GD'] = twosComp2dec(subframe[196:204]) * 2**-31
            eph['IODC'] = int(''.join(subframe[82:84] + subframe[196:204]), 2)
            eph['t_oc'] = int(''.join(subframe[218:234]), 2) * 2**4
            eph['a_f2'] = twosComp2dec(subframe[240:248]) * 2**-55
            eph['a_f1'] = twosComp2dec(subframe[248:264]) * 2**-43
            eph['a_f0'] = twosComp2dec(subframe[270:292]) * 2**-31
            eph['idValid'][0] = 1
        elif subframeID == 2:
            # It contains first part of ephemeris parameters
            eph['IODE_sf2'] = int(''.join(subframe[60:68]), 2)
            eph['C_rs'] = twosComp2dec(subframe[68:84]) * 2**-5
            eph['deltan'] = twosComp2dec(subframe[90:106]) * 2**-43 * gpsPi
            eph['M_0'] = twosComp2dec(subframe[106:114] + subframe[120:144]) * 2**-31 * gpsPi
            eph['C_uc'] = twosComp2dec(subframe[150:166]) * 2**-29
            eph['e'] = int(''.join(subframe[166:174] + subframe[180:204]), 2) * 2**-33
            eph['C_us'] = twosComp2dec(subframe[210:226]) * 2**-29
            eph['sqrtA'] = int(''.join(subframe[226:234] + subframe[240:264]), 2) * 2**-19
            eph['t_oe'] = int(''.join(subframe[270:286]), 2) * 2**4
            eph['idValid'][1] = 2
        elif subframeID == 3:
            # It contains second part of ephemeris parameters
            eph['C_ic'] = twosComp2dec(subframe[60:76]) * 2**-29
            eph['omega_0'] = twosComp2dec(subframe[76:84] + subframe[90:114]) * 2**-31 * gpsPi
            eph['C_is'] = twosComp2dec(subframe[120:136]) * 2**-29
            eph['i_0'] = twosComp2dec(subframe[136:144] + subframe[150:174]) * 2**-31 * gpsPi
            eph['C_rc'] = twosComp2dec(subframe[180:196]) * 2**-5
            eph['omega'] = twosComp2dec(subframe[196:204] + subframe[210:234]) * 2**-31 * gpsPi
            eph['omegaDot'] = twosComp2dec(subframe[240:264]) * 2**-43 * gpsPi
            eph['IODE_sf3'] = int(''.join(subframe[270:278]), 2)
            eph['iDot'] = twosComp2dec(subframe[278:292]) * 2**-43 * gpsPi
            eph['idValid'][2] = 3
        # subframes 4 and 5 not decoded (Almanac, ionospheric model, UTC parameters)

    # Compute the time of week (TOW) of the first sub-frame in the array ====
    # The transmitted TOW is actual TOW of the next subframe; need TOW of the first subframe.
    # Duration of 5 subframes is 30 s, so start time of the first subframe:
    TOW = int(''.join(subframe[30:47]), 2) * 6 - 30
    eph['TOW'] = TOW
    return eph, TOW
