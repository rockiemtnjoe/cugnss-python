import numpy as np

def satpos(transmitTime, prnList, eph):
    """
    SATPOS Calculation of X,Y,Z satellites coordinates at TRANSMITTIME for
    given ephemeris EPH. Coordinates are calculated for each satellite in the
    list PRNLIST.

    [satPositions, satClkCorr] = satpos(transmitTime, prnList, eph)

    Inputs:
        transmitTime  - transmission time: 1 by settings.numberOfChannels
        prnList       - list of PRN-s to be processed
        eph           - ephemeridies of satellites

    Outputs:
        satPositions  - positions of satellites (in ECEF system [X; Y; Z])
        satClkCorr    - correction of satellites clocks in s

    --------------------------------------------------------------------------
                                SoftGNSS v3.0
    --------------------------------------------------------------------------
    Based on Kai Borre 04-09-96
    Copyright (c) by Kai Borre
    Updated by Darius Plausinaitis, Peter Rinder and Nicolaj Bertelsen
    """

    # %% Initialize constants ===================================================
    gpsPi = 3.1415926535898  # Pi used in the GPS coordinate system
    Omegae_dot = 7.2921151467e-5  # Earth rotation rate, [rad/s]
    GM = 3.986005e14  # Earth's universal gravitational constant, [m^3/s^2]
    F = -4.442807633e-10  # Constant, [sec/(meter)^(1/2)]

    numOfSatellites = len(prnList)

    # %% Initialize results =====================================================
    satClkCorr = np.zeros(numOfSatellites)  # Correction of satellites clocks in s
    satPositions = np.zeros((3, numOfSatellites))  # Positions of satellites (ECEF [X; Y; Z])

    # %% Process each satellite =================================================
    for satNr in range(numOfSatellites):
        prn = prnList[satNr]

        # Defensive: skip missing or empty PRNs
        if prn >= len(eph) or not eph[prn] or 't_oc' not in eph[prn]:
            continue

        eph_data = eph[prn-1]

        # %% Find initial satellite clock correction -----------------------------
        #--- Find time difference ---------------------------------------------
        dt = check_t(transmitTime[satNr] - eph_data['t_oc'])

        #--- Calculate clock correction ---------------------------------------
        satClkCorr[satNr] = (eph_data['a_f2'] * dt + eph_data['a_f1']) * dt + \
                             eph_data['a_f0'] - eph_data['T_GD']

        time = transmitTime[satNr] - satClkCorr[satNr]

        # %% Find satellite's position ------------------------------------------
        # Restore semi-major axis
        a = eph_data['sqrtA'] ** 2

        # Time correction
        tk = check_t(time - eph_data['t_oe'])

        # Initial mean motion
        n0 = np.sqrt(GM / a ** 3)
        # Mean motion
        n = n0 + eph_data['deltan']

        # Mean anomaly
        M = eph_data['M_0'] + n * tk
        # Reduce mean anomaly to between 0 and 360 deg
        M = np.remainder(M + 2 * gpsPi, 2 * gpsPi)

        # Initial guess of eccentric anomaly
        E = M

        #--- Iteratively compute eccentric anomaly ----------------------------
        for _ in range(10):
            E_old = E
            E = M + eph_data['e'] * np.sin(E)
            dE = np.remainder(E - E_old, 2 * gpsPi)
            if abs(dE) < 1e-12:
                # Necessary precision is reached, exit from the loop
                break
        # Reduce eccentric anomaly to between 0 and 360 deg
        E = np.remainder(E + 2 * gpsPi, 2 * gpsPi)

        # Relativistic correction
        dtr = F * eph_data['e'] * eph_data['sqrtA'] * np.sin(E)

        # Calculate the true anomaly
        nu = np.arctan2(np.sqrt(1 - eph_data['e'] ** 2) * np.sin(E), np.cos(E) - eph_data['e'])

        # Compute angle phi
        phi = nu + eph_data['omega']
        # Reduce phi to between 0 and 360 deg
        phi = np.remainder(phi, 2 * gpsPi)

        # Correct argument of latitude
        u = phi + eph_data['C_uc'] * np.cos(2 * phi) + eph_data['C_us'] * np.sin(2 * phi)
        # Correct radius
        r = a * (1 - eph_data['e'] * np.cos(E)) + eph_data['C_rc'] * np.cos(2 * phi) + eph_data['C_rs'] * np.sin(2 * phi)
        # Correct inclination
        i = eph_data['i_0'] + eph_data['iDot'] * tk + eph_data['C_ic'] * np.cos(2 * phi) + eph_data['C_is'] * np.sin(2 * phi)

        # SV position in orbital plane
        xk1 = np.cos(u) * r
        yk1 = np.sin(u) * r

        # Compute the angle between the ascending node and the Greenwich meridian
        Omega = eph_data['omega_0'] + (eph_data['omegaDot'] - Omegae_dot) * tk - Omegae_dot * eph_data['t_oe']
        # Reduce to between 0 and 360 deg
        Omega = matlab_rem(Omega + 2 * gpsPi, 2 * gpsPi)

        #--- Compute satellite coordinates ------------------------------------
        xk = xk1 * np.cos(Omega) - yk1 * np.cos(i) * np.sin(Omega)
        yk = xk1 * np.sin(Omega) + yk1 * np.cos(i) * np.cos(Omega)
        zk = yk1 * np.sin(i)

        satPositions[0, satNr] = xk
        satPositions[1, satNr] = yk
        satPositions[2, satNr] = zk

        # %% Include relativistic correction in clock correction ----------------
        satClkCorr[satNr] += dtr  # include relativistic correction

        # %% The following is to calculate sv velocity (currently not used in this version)
        # Computation of SV velocity in ECEF -----------------------------------
        # dE = n/(1-eph_data['e'] * np.cos(E))
        # dphi = np.sqrt(1 - eph_data['e']**2) * dE / (1-eph_data['e'] * np.cos(E))
        # du = dphi + 2*dphi*(-eph_data['C_uc'] * np.sin(2*phi) + eph_data['C_us'] * np.cos(2*phi))
        # dr = a * eph_data['e'] * dE * np.sin(E) + 2*dphi*(-eph_data['C_rc'] * np.sin(2*phi) + eph_data['C_rs'] * np.cos(2*phi))
        # di = eph_data['iDot'] + 2*dphi*(-eph_data['C_ic'] * np.sin(2*phi) + eph_data['C_is'] * np.cos(2*phi))
        # dOmega = eph_data['omegaDot'] - Omegae_dot
        # dxk1 = dr*np.cos(u) - r*du*np.sin(u)
        # dyk1 = dr*np.sin(u) + r*du*np.cos(u)
        # satVolocity[0, satNr] = -yk*dOmega - (dyk1*np.cos(i) - zk*di) * np.sin(Omega) + dxk1*np.cos(Omega)
        # satVolocity[1, satNr] = xk*dOmega  + (dyk1*np.cos(i) - zk*di) * np.cos(Omega) + dxk1*np.sin(Omega)
        # satVolocity[2, satNr] = dyk1*np.sin(i) + yk1*di*np.cos(i)
        # dtrRat = F * eph_data['e'] * eph_data['sqrtA'] * np.cos(E) * dE
        # satClkCorrRat[satNr] = 2* eph_data['a_f2'] * dt + eph_data['a_f1'] + dtrRat

    return satPositions, satClkCorr

def check_t(time):
    """
    Adjust time to be within half a GPS week.
    """
    half_week = 302400  # seconds in half a GPS week
    if time > half_week:
        return time - 2 * half_week
    elif time < -half_week:
        return time + 2 * half_week
    return time

def matlab_rem(a, b):
    """
    MATLAB-style remainder function.
    """
    return a - b * np.fix(a / b)
