import numpy as np

def togeod(a, finv, X, Y, Z):
    """
    TOGEOD   Subroutine to calculate geodetic coordinates latitude, longitude,
             height given Cartesian coordinates X,Y,Z, and reference ellipsoid
             values semi-major axis (a) and the inverse of flattening (finv).

    [dphi, dlambda, h] = togeod(a, finv, X, Y, Z)

      The units of linear parameters X,Y,Z,a must all agree (m,km,mi,ft,..etc)
      The output units of angular quantities will be in decimal degrees
      (15.5 degrees not 15 deg 30 min). The output units of h will be the
      same as the units of X,Y,Z,a.

       Inputs:
           a           - semi-major axis of the reference ellipsoid
           finv        - inverse of flattening of the reference ellipsoid
           X,Y,Z       - Cartesian coordinates

       Outputs:
           dphi        - latitude
           dlambda     - longitude
           h           - height above reference ellipsoid

      Copyright (C) 1987 C. Goad, Columbus, Ohio
      Reprinted with permission of author, 1996
      Fortran code translated into MATLAB
      Kai Borre 03-30-96

    ==========================================================================
    """

    h = 0  # initial value of height
    tolsq = 1.e-10  # convergence tolerance
    maxit = 10  # maximum number of iterations

    # compute radians-to-degree factor
    rtd = 180 / np.pi

    # compute square of eccentricity
    if finv < 1.e-20:
        esq = 0
    else:
        esq = (2 - 1 / finv) / finv

    oneesq = 1 - esq

    # first guess
    # P is distance from spin axis
    P = np.sqrt(X ** 2 + Y ** 2)

    # direct calculation of longitude
    if P > 1.e-20:
        dlambda = np.arctan2(Y, X) * rtd
    else:
        dlambda = 0

    if dlambda < 0:
        dlambda += 360

    # r is distance from origin (0,0,0)
    r = np.sqrt(P ** 2 + Z ** 2)

    if r > 1.e-20:
        sinphi = Z / r
    else:
        sinphi = 0

    dphi = np.arcsin(sinphi)

    # initial value of height = distance from origin minus
    # approximate distance from origin to surface of ellipsoid
    if r < 1.e-20:
        return 0, 0, 0

    h = r - a * (1 - sinphi ** 2 / finv)

    # iterate
    for i in range(maxit):
        sinphi = np.sin(dphi)
        cosphi = np.cos(dphi)

        # compute radius of curvature in prime vertical direction
        N_phi = a / np.sqrt(1 - esq * sinphi ** 2)

        # compute residuals in P and Z
        dP = P - (N_phi + h) * cosphi
        dZ = Z - (N_phi * oneesq + h) * sinphi

        # update height and latitude
        h += (sinphi * dZ + cosphi * dP)
        dphi += (cosphi * dZ - sinphi * dP) / (N_phi + h)

        # test for convergence
        if (dP ** 2 + dZ ** 2) < tolsq:
            break

        # Not Converged--Warn user
        if i == maxit - 1:
            print(f" Problem in TOGEOD, did not converge in {i+1:2.0f} iterations")

    dphi = dphi * rtd  # convert latitude to degrees

    return dphi, dlambda, h