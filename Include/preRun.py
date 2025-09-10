import numpy as np

def pre_run(acq_results, settings):
    """
    Initializes tracking channels from acquisition data. The acquired
    signals are sorted according to the signal strength. This function can be
    modified to use other satellite selection algorithms or to introduce
    acquired signal properties offsets for testing purposes.

    Args:
        acq_results (dict): results from acquisition.
        settings: receiver settings
    Returns:
        List[dict]: List of channel dictionaries with initial tracking setup
    --------------------------------------------------------------------------
    #                           SoftGNSS v3.0
    # 
    # Copyright (C) Darius Plausinaitis
    # Written by Darius Plausinaitis
    # Based on Peter Rinder and Nicolaj Bertelsen
    #--------------------------------------------------------------------------
    #This program is free software; you can redistribute it and/or
    #modify it under the terms of the GNU General Public License
    #as published by the Free Software Foundation; either version 2
    #of the License, or (at your option) any later version.
    #
    #This program is distributed in the hope that it will be useful,
    #but WITHOUT ANY WARRANTY; without even the implied warranty of
    #MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    #GNU General Public License for more details.
    #
    #You should have received a copy of the GNU General Public License
    #along with this program; if not, write to the Free Software
    #Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
    #USA.
    #--------------------------------------------------------------------------
    """
    # Initialize default channel template
    channel_template = {
        "PRN": 0,
        "acquiredFreq": 0,
        "codePhase": 0,
        "status": '-',
    }

    # Create list of default channels
    channels = [channel_template.copy() for _ in range(settings.numberOfChannels)]

    # Sort peaks to find strongest signals, keep the peak index information
    sorted_prns = np.argsort(acq_results["peakMetric"])[::-1]

    # Maximum number of initialized channels is number of detected signals, but
    # not more as the number of channels specified in the settings.
    num_to_init = min(settings.numberOfChannels, np.sum(acq_results["carrFreq"] != 0))

    for ii in range(num_to_init):
        prn = sorted_prns[ii]
        channels[ii]["PRN"] = prn + 1  # MATLAB is 1-indexed
        channels[ii]["acquiredFreq"] = acq_results["carrFreq"][prn]
        channels[ii]["codePhase"] = acq_results["codePhase"][prn]
        channels[ii]["status"] = 'T'  # Tracking state

    return channels
