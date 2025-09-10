def eph_structure_init():
    """
    This is in order to make sure variable 'eph' for each SV has a similar 
    structure when only one or even none of the three requisite messages
    is decoded for a given PRN.

    --------------------------------------------------------------------------
                            CU Multi-GNSS SDR  
    (C) Written by Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
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

    CVS record:
    $Id: ephemeris.m,v 1.1.2.7 2017/03/06 11:38:22 dpl Exp $

    Returns:
        dict: Ephemeris structure initialized with empty/default values
    """

    eph = {
        # Flags for message data decoding. 0 indicates decoding fail, 1 is successful decoding.
        # idValid(1:2) for message type 10 and 11, idValid(3:10) for message type 30:37, idValid(11) for others.
        "idValid": [0] * 5,

        # PRN
        "PRN": [],

        #--- It is subframe 1 -------------------------------------
        # It contains WN, SV clock corrections, health and accuracy
        "weekNumber": [],
        "accuracy": [],
        "health": [],
        "T_GD": [],
        "IODC": [],
        "t_oc": [],
        "a_f2": [],
        "a_f1": [],
        "a_f0": [],

        #--- It is subframe 2 -------------------------------------
        # It contains first part of ephemeris parameters
        "IODE_sf2": [],
        "C_rs": [],
        "deltan": [],
        "M_0": [],
        "C_uc": [],
        "e": [],
        "C_us": [],
        "sqrtA": [],
        "t_oe": [],

        #--- It is subframe 3 -------------------------------------
        # It contains second part of ephemeris parameters
        "C_ic": [],
        "omega_0": [],
        "C_is": [],
        "i_0": [],
        "C_rc": [],
        "omega": [],
        "omegaDot": [],
        "IODE_sf3": [],
        "iDot": [],

        #--- It is subframe 4 -------------------------------------
        # Almanac, ionospheric model, UTC parameters.
        # SV health (PRN: 25-32).
        # Not decoded at the moment.

        #--- It is subframe 5 -------------------------------------
        # SV almanac and health (PRN: 1-24).
        # Almanac reference week number and time.
        # Not decoded at the moment.

        # Tow of first decoded subframe
        "TOW": [],
    }
    return eph
