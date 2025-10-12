import numpy as np

def pre_run(acq_results, settings):
    """
    Initializes tracking channels from acquisition data. The acquired
    signals are sorted according to the signal strength. Only the SVs in LEO
    (with PRN 1 to 32) are tracked - WAAS signals are not tracked.
    
    This function can be
    modified to use other satellite selection algorithms or to introduce
    acquired signal properties offsets for testing purposes.

    Args:
        acq_results (dict): results from acquisition.  Each item in the dictionary
                            is an array.
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

    # Sort peaks to find strongest signals, keep the peak index information
    sorted_idx = np.argsort(acq_results["peakMetric"])[::-1]

    # Now sort through the acquired signals to ensure they are actual 
    # orbiting LEO satellites, not WAAS satellites.

    # Ensure PRNs are a NumPy array before advanced indexing
    prns = np.asarray(acq_results["PRN"])
    sorted_prn_channels = prns[sorted_idx]  

    # Keep only PRNs 1..32 
    allowed_prns = np.arange(1, 33)
    mask = np.isin(sorted_prn_channels, allowed_prns)
    filtered_prns = sorted_prn_channels[mask]      # <- apply the mask to get values
    keep_idx = sorted_idx[mask]

    # Now build a new list of the channel tracking information, with 
    # each element in the list describing the signal to be tracked.
 
    rows = keep_idx[0:settings.numberOfChannels]
    channels = []

    for row in rows:
        channels.append({
            "PRN": int(acq_results["PRN"][row]),
            "acquiredFreq": float(acq_results["carrFreq"][row]),
            "codePhase": int(acq_results["codePhase"][row]),
            "status": "T",
        })
        print (f"PRN{int(acq_results["PRN"][row])} will be tracked - "
               f"peakMetric = {acq_results["peakMetric"][row]:.2f}"
        )

    # pad with defaults if fewer rows than channels
    # This really shouldn't be necessary; subsequent processes should be able to take
    # a list of arbitrary length.
    #
    while len(channels) < settings.numberOfChannels:
        channels.append(channel_template.copy())

    return channels
