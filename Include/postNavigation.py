import numpy as np
from Common.leastSquarePos import least_square_pos
from Include.satpos import satpos
from Common.findUtmZone import find_utm_zone
from Common.calculatePseudoranges import calculate_pseudoranges
from Common.cart2 import cart2geo, cart2utm
from Include.NAVdecoding import NAVdecoding
from copy import deepcopy
from Include.eph_structure_init import eph_structure_init

def postNavigation(trackResults, settings):
    """
    Function calculates navigation solutions for the receiver (pseudoranges,
    positions). At the end it converts coordinates from the WGS84 system to
    the UTM, geocentric or any additional coordinate system.

    [navSolutions, eph] = postNavigation(trackResults, settings)

    Inputs:
        trackResults    - results from the tracking function (structure array).
        settings        - receiver settings.
    Outputs:
        navSolutions    - contains measured pseudoranges, receiver
                        clock error, receiver coordinates in several
                        coordinate systems (at least ECEF and UTM).
        eph             - received ephemerides of all SV (structure array).

    --------------------------------------------------------------------------
                             CU Multi-GNSS SDR  
    (C) Updated by Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
    Based on the original work by Darius Plausinaitis,Peter Rinder, 
    Nicolaj Bertelsen and Dennis M. Akos
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
    $Id: postNavigation.m,v 1.1.2.22 2006/08/09 17:20:11 dpl Exp $
    """

    #--------------------------------------------------------------------------
    # Check is there enough data to obtain any navigation solution ===========
    # It is necessary to have at least three subframes (number 1, 2 and 3) to
    # find satellite coordinates. Then receiver position can be found too.
    # The function requires all 5 subframes, because the tracking starts at
    # arbitrary point. Therefore the first received subframes can be any three
    # from the 5.
    # One subframe length is 6 seconds, therefore we need at least 30 sec long
    # record (5 * 6 = 30 sec = 30000ms). We add extra seconds for the cases,
    # when tracking has started in a middle of a subframe.
    if settings.msToProcess < 36000:
        print('Record is too short. Exiting!')
        return None, None

    #--------------------------------------------------------------------------
    # Pre-allocate space =======================================================
    # Starting positions of the first message in the input bit stream 
    # trackResults.I_P in each channel. The position is PRN code count
    # since start of tracking. Corresponding value will be set to inf 
    # if no valid preambles were detected in the channel.
    subFrameStart = np.full(settings.numberOfChannels, np.inf)

    # Time Of Week (TOW) of the first message(in seconds). Corresponding value
    # will be set to inf if no valid preambles were detected in the channel.
    TOW = np.full(settings.numberOfChannels, np.inf)

    #--- Make a list of channels excluding not tracking channels ---------------
    activeChnList = [i for i, tr in enumerate(trackResults) if tr['status'] != '-']

    #--------------------------------------------------------------------------
    # Decode ephemerides =======================================================
    eph = [deepcopy(eph_structure_init()) for _ in range(33)]

    for ch in activeChnList[:]:
        # Get PRN of current channel
        PRN = trackResults[ch]['PRN']
        print(f'Decoding NAV for PRN {PRN:02d} --------------------')
        # Decode ephemerides and TOW of the first sub-frame
        decoded, subFrameStart[ch], TOW[ch] = NAVdecoding(trackResults[ch]['I_P'], settings)
        eph[PRN - 1].update(decoded)

        # Exclude satellite if it does not have the necessary nav data
        if not (decoded.get('IODC') and decoded.get('IODE_sf2') and decoded.get('IODE_sf3') and decoded.get('health', 0) == 0):
            print(f'    Ephemeris decoding fails for PRN {PRN:02d} !')
            activeChnList.remove(ch)
        else:
            print(f'    Three requisite messages for PRN {PRN:02d} all decoded!')

    #--------------------------------------------------------------------------
    # Check if the number of satellites is still above 3 =====================
    if len(activeChnList) < 4:
        print('Too few satellites with ephemeris data for position calculations. Exiting!')
        return None, None

    #--------------------------------------------------------------------------
    # Set measurement-time point and step  =====================================
    # Find start and end of measurement point locations in IF signal stream with available measurements
    sampleStart = np.array([trackResults[ch]['absoluteSample'][int(subFrameStart[ch])] for ch in activeChnList])
    sampleEnd = np.array([trackResults[ch]['absoluteSample'][-1] for ch in activeChnList])
    # Second term is to make space to avoid index exceeds matrix dimensions, thus a margin of 1 is added.
    sampleStart = float(np.max(sampleStart)) + 1
    sampleEnd = float(np.min(sampleEnd)) - 1

    # Measurement step in unit of IF samples
    step = int(settings.samplingFreq * settings.navSolPeriod / 1000)
    # Number of measurement point from measurement start to end
    measNrSum = int((sampleEnd - sampleStart) / step)

    #--------------------------------------------------------------------------
    # Initialization =========================================================
    # Set the satellite elevations array to INF to include all satellites for
    # the first calculation of receiver position. There is no reference point
    # to find the elevation angle as there is no receiver position estimate at
    # this point.
    satElev = np.full(settings.numberOfChannels, np.inf)

    # Save the active channel list. The list contains satellites that are
    # tracked and have the required ephemeris data. In the next step the list
    # will depend on each satellite's elevation angle, which will change over
    # time.  
    readyChnList = activeChnList[:]

    # Set local time to inf for first calculation of receiver position. After
    # first fix, localTime will be updated by measurement sample step.
    localTime = np.inf

    #--------------------------------------------------------------------------
    #   Do the satellite and receiver position calculations                  #
    #--------------------------------------------------------------------------
    print('Positions are being computed. Please wait...')
    navSolutions = []
    for i in range(measNrSum):
        print(f'Fix: Processing {i+1:02d} of {measNrSum:02d}')
        # Exclude satellites, that are below elevation mask 
        activeNow = [ch for ch in readyChnList if satElev[ch-1] >= settings.elevationMask]

        currSample = sampleStart + step * i
        # Find pseudoranges
        rawP, transmitTime, localTime = calculate_pseudoranges(
            trackResults, subFrameStart, TOW, currSample, localTime, activeNow, settings
        )
        rawP = np.array(rawP)
        transmitTime = np.array(transmitTime)

        # Find satellites positions and clocks corrections
        satPos, satCorr = satpos(transmitTime[activeNow], [trackResults[ch]['PRN'] for ch in activeNow], eph)

        # Find receiver position
        if len(activeNow) > 3:
            # Correct pseudorange for SV clock error
            clkCorr = rawP[activeNow] + satCorr * settings.c
            xyzdt, el, az, dop = least_square_pos(satPos, clkCorr, settings)
            lat, lon, h = cart2geo(*xyzdt[:3], 4)
            utmZone = find_utm_zone(lat, lon)
            e, n, u = cart2utm(*xyzdt[:3], utmZone)

            navSolutions.append({
                'PRN': [trackResults[ch]['PRN'] for ch in activeNow],
                'el': el,
                'az': az,
                'transmitTime': transmitTime[activeNow],
                'satClkCorr': satCorr,
                'rawP': rawP[activeNow],
                'correctedP': rawP[activeNow] + satCorr * settings.c - xyzdt[3],
                'X': xyzdt[0], 'Y': xyzdt[1], 'Z': xyzdt[2], 'dt': xyzdt[3],
                'DOP': dop,
                'latitude': lat, 'longitude': lon, 'height': h,
                'utmZone': utmZone, 'E': e, 'N': n, 'U': u,
                'localTime': localTime - xyzdt[3] / settings.c,
                'currMeasSample': currSample
            })

            # Update the satellites elevations vector
            satElev = el
            # Update local time by measurement step
            localTime += step / settings.samplingFreq
        else:
            print(f'   Measurement No. {i+1}: Not enough information for position solution.')
            # TODO: Known issue. Satellite positions are not updated if the
            # satellites are excluded due to elevation mask. Therefore raising
            # satellites will not be included even if they will be above
            # elevation mask at some point. This would be a good place to
            # update positions of the excluded satellites.
            print('   Exit Program')
            return navSolutions, eph

    return navSolutions, eph
