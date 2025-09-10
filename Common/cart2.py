import numpy as np

def cart2geo(X, Y, Z, i):
    """
    Conversion of Cartesian coordinates (X,Y,Z) to geographical
    coordinates (phi, lambda, h) on a selected reference ellipsoid.

    Choices i of Reference Ellipsoid for Geographical Coordinates:
        0. International Ellipsoid 1924
        1. International Ellipsoid 1967
        2. World Geodetic System 1972
        3. Geodetic Reference System 1980
        4. World Geodetic System 1984

    Returns:
        phi_deg: latitude in degrees
        lambda_deg: longitude in degrees
        h: height above ellipsoid in meters
    """

    # Ellipsoid definitions (semi-major axis and flattening)
    a_list = [6378388, 6378160, 6378135, 6378137, 6378137]
    f_list = [1/297, 1/298.247, 1/298.26, 1/298.257222101, 1/298.257223563]

    a = a_list[i]
    f = f_list[i]

    # Longitude (radians)
    lambda_rad = np.arctan2(Y, X)

    # Intermediate eccentricity squared
    ex2 = (2 - f) * f / ((1 - f) ** 2)
    c = a * np.sqrt(1 + ex2)

    # Initial latitude estimate (radians)
    phi = np.arctan(Z / (np.sqrt(X**2 + Y**2) * (1 - (2 - f) * f)))

    # Iterative computation of latitude and height
    h = 0.1
    oldh = 0
    iterations = 0

    while abs(h - oldh) > 1e-12:
        oldh = h
        N = c / np.sqrt(1 + ex2 * np.cos(phi)**2)
        phi = np.arctan(Z / (np.sqrt(X**2 + Y**2) * (1 - (2 - f) * f * N / (N + h))))
        h = np.sqrt(X**2 + Y**2) / np.cos(phi) - N

        iterations += 1
        if iterations > 100:
            print(f"Warning: did not converge. h-oldh: {h - oldh:.2e}")
            break

    # Convert latitude and longitude to degrees
    phi_deg = np.degrees(phi)
    lambda_deg = np.degrees(lambda_rad)

    return phi_deg, lambda_deg, h

def cart2utm(X, Y, Z, zone):
    """
    Transformation of (X,Y,Z) to (E,N,U) in UTM, zone 'zone'.

    Inputs:
        X, Y, Z - Cartesian coordinates (ITRF96 reference frame)
        zone    - UTM zone of the given position

    Outputs:
        E, N, U - UTM coordinates (Easting, Northing, Uping)
    """

    # Convert zone to integer if it's a string
    if isinstance(zone, str):
        zone = int(''.join(filter(str.isdigit, zone)))

    # Ellipsoid parameters (International Ellipsoid 1924)
    a = 6378388
    f = 1 / 297
    ex2 = (2 - f) * f / ((1 - f) ** 2)
    c = a * np.sqrt(1 + ex2)

    # Coordinate transformation parameters
    alpha = 0.756e-6
    trans = np.array([89.5, 93.8, 127.6])
    scale = 0.9999988

    # Apply transformation to input coordinates
    vec = np.array([X, Y, Z - 4.5])
    R = np.array([[1, -alpha, 0],
                  [alpha, 1, 0],
                  [0, 0, 1]])
    v = scale * R @ vec + trans  # coordinate vector in ED50

    # Initial latitude and longitude estimate (radians)
    L = np.arctan2(v[1], v[0])
    N1 = 6395000  # preliminary value
    B = np.arctan2(v[2] / ((1 - f) ** 2 * N1), np.linalg.norm(v[:2]) / N1)

    # Iterative computation of U
    U = 0.1
    oldU = 0
    while abs(U - oldU) > 1e-4:
        oldU = U
        N1 = c / np.sqrt(1 + ex2 * (np.cos(B) ** 2))
        B = np.arctan2(v[2] / ((1 - f) ** 2 * N1 + U), np.linalg.norm(v[:2]) / (N1 + U))
        U = np.linalg.norm(v[:2]) / np.cos(B) - N1

    # Normalized meridian quadrant (Koenig & Weise)
    m0 = 0.0004
    n = f / (2 - f)
    m = n ** 2 * (1 / 4 + n ** 2 / 64)
    w = (a * (-n - m0 + m * (1 - m0))) / (1 + n)
    Q_n = a + w

    # Easting and longitude of central meridian
    E0 = 500000
    L0 = (zone - 30) * 6 - 3
    L0r = np.deg2rad(L0)

    # Coefficients of trigonometric series (see Koenig & Weise)
    bg = np.array([-3.37077907e-3,
                    4.73444769e-6,
                   -8.29914570e-9,
                    1.58785330e-11])
    gb = np.array([3.37077588e-3,
                   6.62769080e-6,
                   1.78718601e-8,
                   5.49266312e-11])
    gtu = np.array([8.41275991e-4,
                    7.67306686e-7,
                    1.21291230e-9,
                    2.48508228e-12])
    utg = np.array([-8.41276339e-4,
                    -5.95619298e-8,
                    -1.69485209e-10,
                    -2.20473896e-13])

    # Ellipsoidal latitude, longitude to spherical latitude, longitude
    neg_geo = B < 0
    Bg_r = abs(B)
    Bg_r += clsin(bg, 4, 2 * Bg_r)
    Lg_r = L - L0r

    # Spherical latitude, longitude to complementary spherical latitude (N, E)
    cos_BN = np.cos(Bg_r)
    Np = np.arctan2(np.sin(Bg_r), np.cos(Lg_r) * cos_BN)
    Ep = np.arctanh(np.sin(Lg_r) * cos_BN)

    # Spherical normalized N, E to ellipsoidal N, E
    Np *= 2
    Ep *= 2
    dN, dE = clksin(gtu, 4, Np, Ep)
    Np = Np / 2 + dN
    Ep = Ep / 2 + dE

    # Final UTM coordinates
    N = Q_n * Np
    E = Q_n * Ep + E0

    # Southern hemisphere adjustment
    if neg_geo:
        N = -N + 20000000

    return E, N, U

def clsin(coeffs, order, arg):
    """
    Trigonometric series for ellipsoidal to spherical latitude conversion.
    """
    return sum(coeffs[i] * np.sin((i + 1) * arg) for i in range(order))

def clksin(coeffs, order, N, E):
    """
    Trigonometric series for spherical to ellipsoidal N, E conversion.
    """
    dN = sum(coeffs[i] * np.sin((i + 1) * N) * np.cosh((i + 1) * E) for i in range(order))
    dE = sum(coeffs[i] * np.cos((i + 1) * N) * np.sinh((i + 1) * E) for i in range(order))
    return dN, dE
