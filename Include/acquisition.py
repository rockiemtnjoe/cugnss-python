import numpy as np
from Include.makeCaTable import make_ca_table  
from Include.generateCAcode import generate_ca_code  
from Include.plotAcqSearch import plotAcqSearch

def acquisition(long_signal, settings):
    """
    Function performs cold start acquisition on the collected "data". It
    searches for GPS signals of all satellites, which are listed in field
    "acqSatelliteList" in the settings structure. Function saves code phase
    and frequency of the detected signals in the "acqResults" structure.

    Inputs:
        long_signal    - 11 ms of raw signal from the front-end
        settings       - Receiver settings. Provides information about
                         sampling and intermediate frequencies and other
                         parameters including the list of the satellites to
                         be acquired.
    Outputs:
        acqResults     - Function saves code phases and frequencies of the
                         detected signals in the "acqResults" structure. The
                         field "carrFreq" is set to 0 if the signal is not
                         detected for the given PRN number.

    --------------------------------------------------------------------------
                             CU Multi-GNSS SDR
     (C) Updated by Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
     Based on the original work by Darius Plausinaitis,Peter Rinder,
     Nicolaj Bertelsen and Dennis M. Akos
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

    #--------------------------------------------------------------------------
    # Condition input signal to speed up acquisition
    # If input IF signal freq. is too high, a resampling strategy is applied
    # to speed up the acquisition, which is selectable.
    # (Not implemented here, placeholder for now)
    #--------------------------------------------------------------------------

    #--------------------------------------------------------------------------
    # Initialization
    #--------------------------------------------------------------------------
    # Find number of samples per spreading code
    samples_per_code = int(round(settings.samplingFreq / (settings.codeFreqBasis / settings.codeLength)))
    # Find sampling period
    ts = 1 / settings.samplingFreq
    # Find phase points of 2ms local carrier wave (1ms for local duplicate,
    # the other 1ms for zero padding)
    phase_points = np.arange(samples_per_code * 2) * 2 * np.pi * ts
    # Number of the frequency bins for the specified search band
    number_of_freq_bins = int(round(settings.acqSearchBand * 2 / settings.acqSearchStep) + 1)
    # Carrier frequency bins to be searched
    coarse_freq_bin = np.zeros(number_of_freq_bins)

    #--------------------------------------------------------------------------
    # Initialize acqResults
    #--------------------------------------------------------------------------
    # Carrier frequencies of detected signals
    # C/A code phases of detected signals
    # Correlation peak ratios of the detected signals
    numSats = len(settings.acqSatelliteList) + 1
    acqResults = {
        'PRN': np.zeros(numSats),
        'carrFreq': np.zeros(numSats),
        'codePhase': np.zeros(numSats, dtype=int),
        'peakMetric': np.zeros(numSats)
    }
    #--------------------------------------------------------------------------
    # Variables for fine acquisition
    #--------------------------------------------------------------------------
    # Carrier frequency search step for fine acquisition
    fine_search_step = 25
    # Number of the frequency bins for fine acquisition
    num_of_fine_bins = int(round(settings.acqSearchStep / fine_search_step) + 1)
    # Carrier frequencies of the fine frequency bins
    fine_freq_bins = np.zeros(num_of_fine_bins)
    # Phase points of the local carrier wave
    fine_phase_points = np.arange(40 * samples_per_code) * 2 * np.pi * ts

    #--------------------------------------------------------------------------
    # Input signal power for GLRT statistic calculation
    #--------------------------------------------------------------------------

    samples_per_code_val = np.var(long_signal[:samples_per_code], ddof=1)
    sig_power = np.sqrt(samples_per_code_val * samples_per_code)
    sig_power = sig_power

    #--------------------------------------------------------------------------
    # Perform search for all listed PRN numbers ...
    #--------------------------------------------------------------------------

    # TBD: It would be more efficient to compute the frequency bins
    # in advance and store them as an array

    start_freq = settings.IF + settings.acqSearchBand
    stop_freq = settings.IF - settings.acqSearchBand
    num_freq = 2 * int(settings.acqSearchBand / settings.acqSearchStep) + 1
    
    # Create a vector with all the frequencies we will search through
    coarse_freq_bins = np.linspace(start_freq, stop_freq, num_freq)         


    print('(', end='', flush=True)
    for n in range(numSats - 1):
        #--------------------------------------------------------------------------
        # Coarse acquisition
        #--------------------------------------------------------------------------
        # Generate C/A codes and sample them according to the sampling freq.
        ca_codes_table = make_ca_table(settings.acqSatelliteList[n], settings)
        ca_codes_2ms = np.concatenate([ca_codes_table, np.zeros(samples_per_code)])
        # Search results of all frequency bins and code shifts (for one satellite)
        results = np.zeros((number_of_freq_bins, samples_per_code * 2))
        # Perform DFT of C/A code
        ca_code_freq_dom = np.conj(np.fft.fft(ca_codes_2ms))

        # Test the correlation for each frequency bin
        # this is the pythonic way of iterating over all possible frequencies
        # freq_idx goes from 0 to num_freq-1
        # coarse_freq counts through the test frequencies
        for freq_idx, coarse_freq in enumerate(coarse_freq_bins):
            # Generate carrier wave frequency grid
            # coarse_freq_bin[freq_bin_index] = settings.IF + settings.acqSearchBand - settings.acqSearchStep * freq_bin_index
            # Generate local sine and cosine
            sig_carr = np.exp(-1j * coarse_freq * phase_points)
            # Do correlation
            for non_coh_index in range(settings.acqNonCohTime):
                idx_start = non_coh_index * samples_per_code
                idx_end = (non_coh_index + 2) * samples_per_code
                signal = long_signal[idx_start:idx_end]
                if len(signal) < samples_per_code * 2:
                    continue
                # "Remove carrier" from the signal
                I = np.real(sig_carr * signal)
                Q = np.imag(sig_carr * signal)
                # Convert the baseband signal to frequency domain
                IQ_freq_dom = np.fft.fft(I + 1j * Q)
                # Multiplication in the frequency domain (correlation in time domain)
                conv_code_IQ = IQ_freq_dom * ca_code_freq_dom
                # Perform inverse DFT and store correlation results
                coh_result = np.abs(np.fft.ifft(conv_code_IQ))
                # Non-coherent integration
                results[freq_idx, :] += coh_result

        #--------------------------------------------------------------------------
        # Look for correlation peaks for coarse acquisition
        #--------------------------------------------------------------------------
        # Find the correlation peak and the carrier frequency
        coarse_freq_idx = int(np.argmax(np.max(results, axis=1)))
        coarse_freq = coarse_freq_bins[coarse_freq_idx]
        # Find code phase of the same correlation peak
        peak_size = np.max(results)
        code_phase = np.argmax(np.max(results, axis=0))
        # Store GLRT statistic
        acqResults['peakMetric'][n] = peak_size / sig_power / settings.acqNonCohTime

        #
        # Plot the coarse search results
        #
        plotAcqSearch(settings.acqSatelliteList[n],settings, results)
        #--------------------------------------------------------------------------
        # If the result is above threshold, then there is a signal ...
        # Fine carrier frequency search
        #--------------------------------------------------------------------------
        if acqResults['peakMetric'][n] > settings.acqThreshold:
            # Indicate PRN number of the detected signal
            print(f'{settings.acqSatelliteList[n]:02d} ', end='', flush=True)
            # Prepare 20ms code, carrier and input signals
            ca_code = generate_ca_code(settings.acqSatelliteList[n])
            code_value_index = np.floor(ts * np.arange(40 * samples_per_code) / (1 / settings.codeFreqBasis)).astype(int)
            ca_code_40ms = ca_code[np.mod(code_value_index, settings.codeLength)]
            # Take 40ms incoming signal for fine acquisition
            sig_40ms = long_signal[code_phase:code_phase + 40 * samples_per_code]
            # Search different fine freq bins
            fine_result = np.zeros(num_of_fine_bins)
            for fine_bin_index in range(num_of_fine_bins):
                # Carrier frequencies of the frequency bins
                fine_freq_bins[fine_bin_index] = coarse_freq - fine_search_step * fine_bin_index
                # Local carrier signal
                sig_carr_40ms = np.exp(-1j * fine_freq_bins[fine_bin_index] * fine_phase_points)
                # Wipe off code and carrier from incoming signals
                baseband_sig = sig_40ms * ca_code_40ms * sig_carr_40ms
                # Coherent integration for each code
                sum_per_code = np.array([
                    np.sum(baseband_sig[i * samples_per_code:(i + 1) * samples_per_code])
                    for i in range(40)
                ])
                # Search Nav bit edge for 20 cases of Nav bit edge
                max_power = 0
                for com_index in range(20):
                    # Power for 20ms coherent integration
                    com_power = np.abs(np.sum(sum_per_code[com_index:com_index + 20]))
                    # Maximal integration power
                    max_power = max(max_power, com_power)
                fine_result[fine_bin_index] = max_power
            # Find the fine carrier freq.
            max_fin_bin = np.argmax(fine_result)
            acqResults['carrFreq'][n] = fine_freq_bins[max_fin_bin]
            # Save code phase acquisition result
            acqResults['codePhase'][n] = code_phase
            acqResults['PRN'][n] = settings.acqSatelliteList[n]
            # signal found, if IF =0 just change to 1 Hz to allow processing
            if acqResults['carrFreq'][n] == 0:
                acqResults['carrFreq'][n] = 1

        else:
            # No signal with this PRN
            print('. ', end='', flush=True)
    #--------------------------------------------------------------------------
    # Acquisition is over
    #--------------------------------------------------------------------------
    print(')')
    return acqResults
