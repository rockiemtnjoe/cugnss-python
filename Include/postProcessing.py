"""
Script postProcessing.py processes the raw signal from the specified data
file (in settings) operating on blocks of 37 seconds of data.

First it runs acquisition code identifying the satellites in the file,
then the code and carrier for each of the satellites are tracked, storing
the 1msec accumulations. After processing all satellites in the 37 sec
data block, then postNavigation is called. It calculates pseudoranges
and attempts a position solution. At the end plots are made for that
block of data.

--------------------------------------------------------------------------
                        CU Multi-GNSS SDR  
(C) Updated by Jakob Almqvist, Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
Based on the original work by Darius Plausinaitis, Peter Rinder, 
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

                        THE SCRIPT "RECIPE"

The purpose of this script is to combine all parts of the software
receiver.

1.1) Open the data file for the processing and seek to desired point.

2.1) Acquire satellites

3.1) Initialize channels (preRun).
3.2) Pass the channel structure and the file identifier to the tracking
     function. It will read and process the data. The tracking results are
     stored in the trkResults structure. The results can be accessed this
     way (the results are stored each millisecond):
     trkResults[channelNumber].XXX[fromMillisecond : toMillisecond], where
     XXX is a field name of the result (e.g. I_P, codePhase etc.)

4)   Pass tracking results to the navigation solution function. It will
     decode navigation messages, find satellite positions, measure
     pseudoranges and find receiver position.

5)   Plot the results.
"""
import numpy as np
from datetime import datetime
from Include.acquisition import acquisition
from Include.preRun import pre_run
from Include.tracking import tracking
from Include.postNavigation import postNavigation
from Include.showChannelStatus import show_channel_status
from Include.plotAcquisition import plotAcquisition
from Include.plotTracking import plotTracking
from Include.plotNavigation import plotNavigation

def postProcessing(settings):
    # %% Initialization =========================================================
    print('Starting processing...')

    try:
        fid = open(settings.fileName, 'rb')
    except Exception as e:
        # Error while opening the data file.
        raise RuntimeError(f"Unable to read file {settings.fileName}: {e}")

    # Initialize the multiplier to adjust for the data type
    data_adapt_coeff = 1 if settings.fileType == 1 else 2

    # Move the starting point of processing. Can be used to start the
    # signal processing at any point in the data record (e.g. good for long
    # records or for signal processing in blocks).
    fid.seek(data_adapt_coeff * settings.skipNumberOfBytes, 0)

    # %% Acquisition ============================================================
    acqResults = None
    if settings.skipAcquisition == 0:
        # Find number of samples per spreading code
        samples_per_code = int(round(settings.samplingFreq / (settings.codeFreqBasis / settings.codeLength)))
        # At least 42ms of signal are needed for fine frequency estimation
        code_len = max(42, settings.acqNonCohTime + 2)
        num_samples = data_adapt_coeff * code_len * samples_per_code

        # Read data for acquisition.
        if settings.dataType == 'schar':
            dtype = np.int8
        elif settings.dataType == 'short':
            dtype = np.int16
        elif settings.dataType == 'float':
            dtype = np.float32
        else:
            raise ValueError(f"Unsupported dataType: {settings.dataType}")
        data = np.fromfile(fid, dtype=dtype, count=num_samples)
        if data.size < num_samples:
            raise ValueError('Could not read enough data from the data file.')

        if data_adapt_coeff == 2:
            # For complex data, separate I and Q
            data_i = data[::2]
            data_q = data[1::2]
            data = data_i + 1j * data_q

        # --- Do the acquisition -------------------------------------------
        print('   Acquiring satellites...')
        acqResults = acquisition(data, settings)

        # Save acquisition results
        np.save("acqResults.npy", acqResults)
        print('   Acquisition results saved to acqResults.npy')
        
    else:
        # Load the npy file
        try:
            acqResults = np.load("acqResults.npy", allow_pickle=True).item()
            print('   Loaded acquisition results from acqResults.npy')
        except FileNotFoundError:
            print('   acqResults.npy not found. Please run acquisition first.')
            fid.close()
            return
        except Exception as e:
            print(f'   Error loading acqResults.npy: {e}')
            fid.close()
            return


    # %% Initialize channels and prepare for the run ============================
    # Start further processing only if a GNSS signal was acquired (the
    # field FREQUENCY will be set to 0 for all not acquired signals)
    if acqResults is not None and np.any(acqResults['carrFreq']):
        channel = pre_run(acqResults, settings)
        show_channel_status(channel, settings)
    else:
        # No satellites to track, exit
        print('No GNSS signals detected, signal processing finished.')
        fid.close()
        return

    # %% Track the signal =======================================================
    start_time = datetime.now()
    print(f'   Tracking started at {start_time.strftime("%Y-%m-%d %H:%M:%S")}')
    # Process all channels for given data block
    trkResults, _ = tracking(fid, channel, settings)
    fid.close()
    
    elapsed = datetime.now() - start_time
    print(f'   Tracking is over (elapsed time {elapsed})')

    np.save("trkResults.npy", trkResults)
    print('   Tracking results saved to trkResults.npy')

    # %% Calculate navigation solutions =========================================
    print('   Calculating navigation solutions...')
    navResults, eph = postNavigation(trkResults, settings)
    np.save("navResults.npy", navResults)
    print('   Navigation results saved to navResults.npy')

    print('   Processing is complete for this data block')
    print('Post processing of the signal is over.')

    # %% Plot all results ===================================================
    print('   Plotting results...')
    if settings.plotAcquisition:
        plotAcquisition(acqResults)
    if settings.plotTracking:
        plotTracking(range(1, settings.numberOfChannels + 1), trkResults, settings)
    if settings.plotNavigation:
        plotNavigation(navResults, settings)
    print('Post processing of the signal is over.')
