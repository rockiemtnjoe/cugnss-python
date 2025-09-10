import numpy as np
from .e_r_corr import e_r_corr
from Include.topocent import topocent
from .tropo import tropo

def least_square_pos(satpos, obs, settings):
    # Function calculates the Least Square Solution.
    #
    # [pos, el, az, dop] = least_square_pos(satpos, obs, settings)
    #
    #   Inputs:
    #       satpos      - Satellites positions (in ECEF system: [X; Y; Z] -
    #                   one column per satellite)
    #       obs         - Observations - the pseudorange measurements to each
    #                   satellite corrected by SV clock error
    #                   (e.g. [20000000 21000000 .... .... .... .... ....]) 
    #       settings    - receiver settings
    #
    #   Outputs:
    #       pos         - receiver position and receiver clock error 
    #                   (in ECEF system: [X, Y, Z, dt]) 
    #       el          - Satellites elevation angles (degrees)
    #       az          - Satellites azimuth angles (degrees)
    #       dop         - Dilutions Of Precision ([GDOP PDOP HDOP VDOP TDOP])
    # --------------------------------------------------------------------------
    #                           SoftGNSS v3.0
    # --------------------------------------------------------------------------
    # Based on Kai Borre
    # Copyright (c) by Kai Borre
    # Updated by Darius Plausinaitis, Peter Rinder and Nicolaj Bertelsen
    # CVS record:
    # $Id: leastSquarePos.m,v 1.1.2.12 2006/08/22 13:45:59 dpl Exp $
    # ==========================================================================

    # === Initialization =======================================================
    n_iter = 10
    dtr = np.pi / 180
    pos = np.zeros(4)  # center of earth
    X = satpos
    n_sats = X.shape[1]

    A = np.zeros((n_sats, 4))
    omc = np.zeros(n_sats)
    az = np.zeros(n_sats)
    el = np.zeros(n_sats)

    # === Iteratively find receiver position ===================================
    for it in range(n_iter):
        for i in range(n_sats):
            if it == 0:
                # --- Initialize variables at the first iteration --------------
                Rot_X = X[:, i]
                trop = 2.0
            else:
                # --- Update equations -----------------------------------------
                rho2 = np.sum((X[:, i] - pos[:3])**2)
                traveltime = np.sqrt(rho2) / settings.c

                # --- Correct satellite position (due to earth rotation) -------
                # Convert SV position at signal transmitting time to position 
                # at signal receiving time. ECEF always changes with time as 
                # earth rotates.
                Rot_X = e_r_corr(traveltime, X[:, i])

                # --- Find the elevation angle of the satellite ----------------
                az[i], el[i], _ = topocent(pos[:3], Rot_X - pos[:3])

                if getattr(settings, "useTropCorr", 0) == 1:
                    # --- Calculate tropospheric correction --------------------
                    trop = tropo(np.sin(el[i] * dtr), 0.0, 1013.0, 293.0, 50.0, 0.0, 0.0, 0.0)
                else:
                    # Do not calculate or apply the tropospheric corrections
                    trop = 0.0

            # --- Apply the corrections ----------------------------------------
            omc[i] = obs[i] - np.linalg.norm(Rot_X - pos[:3]) - pos[3] - trop

            # --- Construct the A matrix ---------------------------------------
            r = np.linalg.norm(Rot_X - pos[:3])
            if r == 0:
                A[i, :] = 0  # Avoid division by zero
            else:
                A[i, :] = [-(Rot_X[0] - pos[0]) / r,
                           -(Rot_X[1] - pos[1]) / r,
                           -(Rot_X[2] - pos[2]) / r,
                           1]

        # These lines allow the code to exit gracefully in case of any errors
        if np.linalg.matrix_rank(A) < 4:
            print("Cannot get a converged solution!")
            return np.zeros(4), el, az, np.inf * np.ones(5)

        # --- Find position update (in the least squares sense)-----------------
        x = np.linalg.lstsq(A, omc, rcond=None)[0]

        # --- Apply position update --------------------------------------------
        pos += x

    # === Calculate Dilution Of Precision ======================================
    Q = np.linalg.inv(A.T @ A)
    dop = np.zeros(5)
    dop[0] = np.sqrt(np.trace(Q))                  # GDOP    
    dop[1] = np.sqrt(Q[0, 0] + Q[1, 1] + Q[2, 2])   # PDOP
    dop[2] = np.sqrt(Q[0, 0] + Q[1, 1])             # HDOP
    dop[3] = np.sqrt(Q[2, 2])                       # VDOP
    dop[4] = np.sqrt(Q[3, 3])                       # TDOP

    return pos, el, az, dop
