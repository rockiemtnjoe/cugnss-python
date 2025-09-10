import numpy as np

def cno_vsm(I, Q, T):
    """
    Calculate C/No using the Variance Summing Method

    Parameters
    ----------
    I : np.ndarray
        Prompt In Phase values of the signal from Tracking
    Q : np.ndarray
        Prompt Quadrature Phase values of the signal from Tracking
    T : float
        Accumulation interval in Tracking (in seconds)

    Returns
    -------
    float
        Estimated C/No for the given values of I and Q (in dB-Hz)

    Notes
    -----
    This function implements the Variance Summing Method for estimating
    Carrier-to-Noise ratio (C/No) as described in the original MATLAB code.

    Steps:
        1. Calculate Power: Z = I**2 + Q**2
        2. Calculate mean and variance of Power: Zm = mean(Z), Zv = var(Z)
        3. Calculate average carrier power: Pav = sqrt(Zm**2 - Zv)
        4. Calculate variance of the noise: Nv = 0.5 * (Zm - Pav)
        5. Calculate C/No: CNo = 10 * log10(abs((1/T) * Pav / (2*Nv)))

    Copyright (C) D.M.Akos
    Written by Sirish Jetti

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
    """
    # Calculate Power
    Z = I**2 + Q**2
    # Calculate the mean and variance of the Power
    Zm = np.mean(Z)
    Zv = np.var(Z)
    # Calculate the average carrier power
    diff = Zm**2 - Zv
    if diff <= 0:
        return float('nan')
    Pav = np.sqrt(diff)
    # Calculate the variance of the noise
    Nv = 0.5 * (Zm - Pav)
    if Nv <= 0:
        return float('nan')
    # Calculate C/No
    CNo = 10 * np.log10(abs((1 / T) * Pav / (2 * Nv)))
    return CNo
