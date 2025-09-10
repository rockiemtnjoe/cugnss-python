import numpy as np
from Common.togeod import togeod

def topocent(X, dx):
    # TOPOCENT  Transformation of vector dx into topocentric coordinate
    #           system with origin at X.
    #           Both parameters are 3 by 1 vectors.
    #
    # [Az, El, D] = topocent(X, dx)
    #
    #   Inputs:
    #       X           - vector origin coordinates (in ECEF system [X, Y, Z])
    #       dx          - vector ([dX, dY, dZ])
    #
    #   Outputs:
    #       D           - vector length. Units like units of the input
    #       Az          - azimuth from north positive clockwise, degrees
    #       El          - elevation angle, degrees
    #
    # Kai Borre 11-24-96
    # Copyright (c) by Kai Borre
    #
    # CVS record:
    # $Id: topocent.m,v 1.1.1.1.2.4 2006/08/22 13:45:59 dpl Exp $
    # ==========================================================================

    dtr = np.pi / 180

    # Ensure input vectors are 1D arrays of length 3
    X = np.asarray(X)
    if X.ndim == 2 and X.shape[0] == 1:
        X = X.flatten()
    if X.size != 3:
        raise ValueError("X must be a vector of length 3 (ECEF coordinates).")
    dx = np.asarray(dx)
    if dx.ndim == 2 and dx.shape[0] == 1:
        dx = dx.flatten()
    if dx.size != 3:
        raise ValueError("dx must be a vector of length 3 (ECEF delta).")

    # Geodetic coordinates of origin
    phi, lambda_, _ = togeod(6378137, 298.257223563, X[0], X[1], X[2])

    # Trigonometric functions
    cl = np.cos(lambda_ * dtr)
    sl = np.sin(lambda_ * dtr)
    cb = np.cos(phi * dtr)
    sb = np.sin(phi * dtr)

    # Transformation matrix
    F = np.array([[-sl, -sb * cl, cb * cl],
                  [ cl, -sb * sl, cb * sl],
                  [  0,     cb,      sb   ]])

    # Transform dx to local topocentric coordinates
    local_vector = F.T @ dx
    E = local_vector[0]
    N = local_vector[1]
    U = local_vector[2]

    # Horizontal distance
    hor_dis = np.sqrt(E**2 + N**2)

    # Azimuth and elevation calculation
    if hor_dis < 1e-20:
        Az = 0
        El = 90
    else:
        Az = np.degrees(np.arctan2(E, N))
        El = np.degrees(np.arctan2(U, hor_dis))

    # Ensure azimuth is positive
    if Az < 0:
        Az += 360

    # Vector length
    D = np.linalg.norm(dx)
    # %%%%%%%%% end topocent.py %%%%%%%%%
    return Az, El, D