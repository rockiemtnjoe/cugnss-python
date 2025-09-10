import numpy as np

def e_r_corr(travel_time, x_sat):
    #E_R_CORR  Returns rotated satellite ECEF coordinates due to Earth
    #rotation during signal travel time
    #
    #x_sat_rot = e_r_corr(travel_time, x_sat)
    #
    #   Inputs:
    #       travel_time  - signal travel time
    #       x_sat        - satellite's ECEF coordinates
    #
    #   Outputs:
    #       x_sat_rot    - rotated satellite's coordinates (ECEF)
    #
    #Written by Kai Borre
    #Copyright (c) by Kai Borre
    #
    # CVS record:
    # $Id: e_r_corr.m,v 1.1.1.1.2.6 2006/08/22 13:45:59 dpl Exp $
    #==========================================================================

    Omegae_dot = 7.292115147e-5  # rad/sec

    #--- Find rotation angle --------------------------------------------------
    omegatau = Omegae_dot * travel_time

    #--- Make a rotation matrix -----------------------------------------------
    R3 = np.array([
        [np.cos(omegatau),  np.sin(omegatau), 0],
        [-np.sin(omegatau), np.cos(omegatau), 0],
        [0,                 0,                1]
    ])

    #--- Do the rotation ------------------------------------------------------
    x_sat_rot = R3 @ x_sat

    return x_sat_rot

#%%%%%%% end e_r_corr.py %%%%%%%%%%%%%%%%%%%%
