def find_utm_zone(latitude, longitude):
    """
    Function finds the UTM zone number for given longitude and latitude.
    The longitude value must be between -180 (180 degree West) and 180 (180 degree East) degree.
    The latitude must be within -80 (80 degree South) and 84 (84 degree North).

    Usage:
        utm_zone = find_utm_zone(latitude, longitude)

    Latitude and longitude must be in decimal degrees (e.g. 15.5 degrees not 15 deg 30 min).

    --------------------------------------------------------------------------
                                SoftGNSS v3.0

    Copyright (C) Darius Plausinaitis
    Written by Darius Plausinaitis
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
    ==========================================================================

    CVS record:
    $Id: findUtmZone.m,v 1.1.2.2 2006/08/22 13:45:59 dpl Exp $
    """

    # Check value bounds =====================================================
    if longitude > 180 or longitude < -180:
        raise ValueError("Longitude value exceeds limits (-180 to 180).")

    if latitude > 84 or latitude < -80:
        raise ValueError("Latitude value exceeds limits (-80 to 84).")

    # Find zone ==============================================================
    # Start at 180 deg west = -180 deg
    utm_zone = int((180 + longitude) / 6) + 1

    # Correct zone numbers for particular areas ==============================
    if latitude > 72:
        # Corrections for zones 31 33 35 37
        if 0 <= longitude < 9:
            utm_zone = 31
        elif 9 <= longitude < 21:
            utm_zone = 33
        elif 21 <= longitude < 33:
            utm_zone = 35
        elif 33 <= longitude < 42:
            utm_zone = 37

    elif 56 <= latitude < 64:
        # Correction for zone 32
        if 3 <= longitude < 12:
            utm_zone = 32

    return str(utm_zone)