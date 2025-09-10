"""
calculate_pseudoranges finds relative pseudoranges for all satellites
listed in channel_list at the specified millisecond of the processed
signal. The pseudoranges contain unknown receiver clock offset. It can be
found by the least squares position search procedure.

Inputs:
    track_results    - output from the tracking function
    sub_frame_start  - array contains positions of the first preamble in each channel.
                       The position is ms count since start of tracking. Corresponding
                       value will be set to 0 if no valid preambles were detected in
                       the channel. 1 by settings.numberOfChannels
    TOW              - Time Of Week (TOW) of the first sub-frame in the bit stream (in seconds)
    curr_meas_sample - current measurement sample location (measurement time)
    local_time       - local time (in GPST) at measurement time
    channel_list     - list of channels to be processed
    settings         - receiver settings

Outputs:
    pseudoranges     - relative pseudoranges to the satellites
    transmit_time    - transmitting time of channels to be processed corresponding to measurement time
    local_time       - local time (in GPST) at measurement time

--------------------------------------------------------------------------
CU Multi-GNSS SDR
(C) Written by Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
based on the original work by Darius Plausinaitis,
Peter Rinder, Nicolaj Bertelsen and Dennis M. Akos
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
--------------------------------------------------------------------------
"""

def calculate_pseudoranges(track_results, sub_frame_start, TOW, curr_meas_sample, local_time, channel_list, settings):
    """
    Calculate relative pseudoranges and transmit times for satellites.

    Args:
        track_results: List of tracking result objects per channel.
        sub_frame_start: List of start indices of subframes for each channel.
        TOW: List of time-of-week values for each channel.
        curr_meas_sample: Current measurement sample index.
        local_time: Current local receiver time.
        channel_list: List of active channel indices.
        settings: Configuration settings object.

    Returns:
        pseudoranges: List of computed pseudoranges.
        transmit_time: List of satellite transmit times.
        local_time: Updated local time.
    """

    # Transmitting Time of all channels at current measurement sample location
    transmit_time = [float('inf')] * settings.numberOfChannels

    # For all channels in the list ...
    for channel_nr in channel_list:
        # Find index of I_P stream whose integration contains current measurement point location
        index = 0
        for i, sample in enumerate(track_results[channel_nr]['absoluteSample']):
            if sample > curr_meas_sample:
                index = i
                break
        index -= 1

        # Update the phase step based on code freq and sampling frequency
        code_freq = track_results[channel_nr]['codeFreq'][index]
        code_phase_step = code_freq / settings.samplingFreq

        # Code phase from start of a PRN code to current measurement sample location
        rem_code_phase = track_results[channel_nr]['remCodePhase'][index]
        abs_sample = track_results[channel_nr]['absoluteSample'][index]
        code_phase = rem_code_phase + code_phase_step * (curr_meas_sample - abs_sample)

        # Transmitting Time (in unit of s) at current measurement sample location
        # code_phase/settings.codeLength: fraction part of a PRN code
        # index - sub_frame_start[channel_nr]: integer number of PRN code
        transmit_time[channel_nr] = (
            (code_phase / settings.codeLength + index - sub_frame_start[channel_nr]) *
            settings.codeLength / settings.codeFreqBasis + TOW[channel_nr]
        )

    # At first time of fix, local time is initialized by transmit_time and settings.startOffset
    if local_time == float('inf'):
        max_time = max([transmit_time[i] for i in channel_list])
        local_time = max_time + settings.startOffset / 1000.0

    # Convert travel time to a distance
    # The speed of light must be converted from meters per second to meters per millisecond.
    pseudoranges = [(local_time - t) * settings.c for t in transmit_time]

    return pseudoranges, transmit_time, local_time
