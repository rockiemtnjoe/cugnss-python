import numpy as np
from .generateCAcode import generate_ca_code
from Common.CNoVSM import cno_vsm
from Common.calcLoopCoef import calc_loop_coef
import sys

def print_progress(loopCnt, codePeriods, channelNr, PRN, CNo):
    # Progress bar for tracking status
    bar_len = 40
    percent = loopCnt / codePeriods
    filled_len = int(round(bar_len * percent))
    bar = '█' * filled_len + '-' * (bar_len - filled_len)
    sys.stdout.write(
        f'\rTracking: Ch {channelNr + 1} PRN {PRN} |{bar}| {int(percent * 100)}% '
        f'({loopCnt}/{codePeriods} ms) C/No: {CNo if isinstance(CNo, str) else round(CNo,1)} dB-Hz'
    )
    sys.stdout.flush()
    if loopCnt == codePeriods:
        print()

def tracking(fid, channel, settings):
    """
    Performs code and carrier tracking for all channels.

    Inputs:
        fid      - file identifier of the signal record.
        channel  - PRN, carrier frequencies and code phases of all satellites to be tracked.
        settings - receiver settings.

    Outputs:
        trackResults - tracking results (structure array). Contains in-phase prompt outputs and absolute spreading
                       code's starting positions, together with other observation data from the tracking loops.
                       All are saved every millisecond.
        channel      - updated channel information.
    """

    #--------------------------------------------------------------------------
    #                         CU Multi-GNSS SDR  
    # (C) Updated by Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
    # Based on the original work by Darius Plausinaitis,Peter Rinder, 
    # Nicolaj Bertelsen and Dennis M. Akos
    #--------------------------------------------------------------------------

    # This program is free software; you can redistribute it and/or
    # modify it under the terms of the GNU General Public License
    # as published by the Free Software Foundation; either version 2
    # of the License, or (at your option) any later version.
    #
    # This program is distributed in the hope that it will be useful,
    # but WITHOUT ANY WARRANTY; without even the implied warranty of
    # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    # GNU General Public License for more details.
    #
    # You should have received a copy of the GNU General Public License
    # along with this program; if not, write to the Free Software
    # Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
    # USA.
    #--------------------------------------------------------------------------

    # === Initialize result structure ============================================
    def init_result():
        # Channel status: No tracked signal, or lost lock
        return {
            'status': '-',  # No tracked signal, or lost lock
            # The absolute sample in the record of the C/A code start:
            'absoluteSample': np.zeros(settings.msToProcess),
            # Freq of the PRN code:
            'codeFreq': np.full(settings.msToProcess, np.inf),
            # Frequency of the tracked carrier wave:
            'carrFreq': np.full(settings.msToProcess, np.inf),
            # Outputs from the correlators (In-phase):
            'I_P': np.zeros(settings.msToProcess),
            'I_E': np.zeros(settings.msToProcess),
            'I_L': np.zeros(settings.msToProcess),
            # Outputs from the correlators (Quadrature-phase):
            'Q_E': np.zeros(settings.msToProcess),
            'Q_P': np.zeros(settings.msToProcess),
            'Q_L': np.zeros(settings.msToProcess),
            # Loop discriminators
            'dllDiscr': np.full(settings.msToProcess, np.inf),
            'dllDiscrFilt': np.full(settings.msToProcess, np.inf),
            'pllDiscr': np.full(settings.msToProcess, np.inf),
            'pllDiscrFilt': np.full(settings.msToProcess, np.inf),
            # Remain code and carrier phase
            'remCodePhase': np.full(settings.msToProcess, np.inf),
            'remCarrPhase': np.full(settings.msToProcess, np.inf),
            # C/No
            'CNo': {
                'VSMValue': np.zeros(settings.msToProcess // settings.CNo.VSMinterval),
                'VSMIndex': np.zeros(settings.msToProcess // settings.CNo.VSMinterval)
            },
            # PRN number
            'PRN': 0
        }

    # Copy initial settings for all channels
    trackResults = [init_result() for _ in range(settings.numberOfChannels)]

    # === Initialize tracking variables ==========================================
    # Signal period to be processed
    codePeriods = settings.msToProcess  # For GPS one C/A code is one ms

    # --- DLL variables ---------------------------------------------------------
    # Define early-late offset (in chips)
    earlyLateSpc = settings.dllCorrelatorSpacing
    # Summation interval
    PDIcode = settings.intTime
    # Calculate filter coefficient values
    tau1code, tau2code = calc_loop_coef(settings.dllNoiseBandwidth, settings.dllDampingRatio, 1.0)

    # --- PLL variables ---------------------------------------------------------
    # Summation interval
    PDIcarr = settings.intTime
    # Calculate filter coefficient values
    tau1carr, tau2carr = calc_loop_coef(settings.pllNoiseBandwidth, settings.pllDampingRatio, 0.25)

    # Data adaptation coefficient (1 for real, 2 for complex)
    dataAdaptCoeff = 1 if settings.fileType == 1 else 2

    # === Start processing channels ==============================================
    for channelNr in range(settings.numberOfChannels):

        ch = channel[channelNr]
        # Only process if PRN is non zero (acquisition was successful)
        if ch['PRN'] == 0:
            continue

        # Save additional information - each channel's tracked PRN
        trackResults[channelNr]['PRN'] = ch['PRN']

        # Move the starting point of processing. Can be used to start the
        # signal processing at any point in the data record (e.g. for long
        # records). In addition skip through that data file to start at the
        # appropriate sample (corresponding to code phase).
        bytes_per_element = 2 if settings.dataType == 'int16' else 1
        seek_offset = (settings.skipNumberOfBytes + ch['codePhase'] * dataAdaptCoeff) * bytes_per_element
        fid.seek(seek_offset, 0)

        # Get a vector with the C/A code sampled 1x/chip
        caCode = generate_ca_code(ch['PRN'])
        # Then make it possible to do early and late versions
        caCode = np.concatenate([[caCode[-1]], caCode, [caCode[0]]])

        # --- Perform various initializations ------------------------------------
        # Define initial code frequency basis of NCO
        codeFreq = settings.codeFreqBasis
        # Define residual code phase (in chips)
        remCodePhase = 0.0
        # Define carrier frequency which is used over whole tracking period
        carrFreq = ch['acquiredFreq']
        carrFreqBasis = carrFreq
        # Define residual carrier phase
        remCarrPhase = 0.0

        # Code tracking loop parameters
        oldCodeNco = oldCodeError = 0.0
        # Carrier/Costas loop parameters
        oldCarrNco = oldCarrError = 0.0
        # C/No computation
        vsmCnt = 0
        CNo = 0

        # === Process the number of specified code periods =======================
        for loopCnt in range(settings.msToProcess):

            # --- Progress Bar / GUI update --------------------------------------
            # The progress bar is updated every 50ms.
            if loopCnt % 50 == 0 or loopCnt == settings.msToProcess - 1:
                print_progress(loopCnt + 1, settings.msToProcess, channelNr, ch['PRN'], CNo)

            # Record sample number (based on 8bit samples)
            trackResults[channelNr]['absoluteSample'][loopCnt] = fid.tell() / dataAdaptCoeff / (2 if settings.dataType == 'int16' else 1)

            # Update the phasestep based on code freq (variable) and sampling frequency (fixed)
            codePhaseStep = codeFreq / settings.samplingFreq

            # Find the size of a "block" or code period in whole samples
            blksize = int(np.ceil((settings.codeLength - remCodePhase) / codePhaseStep))

            # Read in the appropriate number of samples to process this iteration
            bytes_per_element = 2 if settings.dataType == 'int16' else 1
            num_elements = dataAdaptCoeff * blksize
            num_bytes = num_elements * bytes_per_element

            raw_bytes = fid.read(num_bytes)
            if len(raw_bytes) < num_bytes:
                print(f"Not enough bytes read on channel {channelNr}. Expected {num_bytes}, got {len(raw_bytes)}. Exiting.")
                return trackResults, channel

            raw_array = np.frombuffer(raw_bytes, dtype=np.int16 if settings.dataType == 'int16' else np.int8)

            # For complex data
            if dataAdaptCoeff == 2:
                # check even number of elements
                if len(raw_array) % 2 != 0:
                    raw_array = raw_array[:-1]
                rawSignal = raw_array[::2] + 1j * raw_array[1::2]
            else:
                rawSignal = raw_array.astype(np.float64)

            # --- Set up all the code phase tracking information ------------------
            # Save remCodePhase for current correlation
            trackResults[channelNr]['remCodePhase'][loopCnt] = remCodePhase

            # Define index into early code vector
            def ca_slice(tcode):
                idx = np.ceil(tcode).astype(int)
                return caCode[np.clip(idx, 0, len(caCode) - 1)]

            tcode = remCodePhase - earlyLateSpc + np.arange(blksize) * codePhaseStep
            earlyCode = ca_slice(tcode)

            tcode = remCodePhase + earlyLateSpc + np.arange(blksize) * codePhaseStep
            lateCode = ca_slice(tcode)

            tcode = remCodePhase + np.arange(blksize) * codePhaseStep
            promptCode = ca_slice(tcode)

            # Remaining code phase for each tracking update
            remCodePhase = (tcode[-1] + codePhaseStep) - settings.codeLength

            # --- Generate the carrier frequency to mix the signal to baseband ----
            # Save remCarrPhase for current correlation
            trackResults[channelNr]['remCarrPhase'][loopCnt] = remCarrPhase

            # Get the argument to sin/cos functions
            time = np.arange(blksize + 1) / settings.samplingFreq
            trigarg = 2.0 * np.pi * carrFreq * time + remCarrPhase
            # Remaining carrier phase for each tracking update
            remCarrPhase = np.remainder(trigarg[-1], 2 * np.pi)

            # Finally compute the signal to mix the collected data to baseband
            carrsig = np.exp(-1j * trigarg[:-1])

            # --- Do correlation to Generate the six standard accumulated values ---
            # First mix to baseband
            iBaseband = np.real(carrsig * rawSignal[:blksize])
            qBaseband = np.imag(carrsig * rawSignal[:blksize])

            # Now get early, late, and prompt values for each
            I_E, Q_E = np.sum(earlyCode * iBaseband), np.sum(earlyCode * qBaseband)
            I_P, Q_P = np.sum(promptCode * iBaseband), np.sum(promptCode * qBaseband)
            I_L, Q_L = np.sum(lateCode * iBaseband), np.sum(lateCode * qBaseband)

            # --- Find PLL error and update carrier NCO ---------------------------
            # Implement carrier loop discriminator (phase detector)
            carrError = np.arctan(Q_P / I_P) / (2.0 * np.pi)
            # Implement carrier loop filter and generate NCO command
            carrNco = oldCarrNco + (tau2carr / tau1carr) * (carrError - oldCarrError) + carrError * (PDIcarr / tau1carr)
            oldCarrNco, oldCarrError = carrNco, carrError

            # Save carrier frequency for current correlation
            trackResults[channelNr]['carrFreq'][loopCnt] = carrFreq
            # Modify carrier freq based on NCO command
            carrFreq = carrFreqBasis + carrNco

            # --- Find DLL error and update code NCO ------------------------------
            codeError = (np.hypot(I_E, Q_E) - np.hypot(I_L, Q_L)) / (np.hypot(I_E, Q_E) + np.hypot(I_L, Q_L))
            # Implement code loop filter and generate NCO command
            codeNco = oldCodeNco + (tau2code / tau1code) * (codeError - oldCodeError) + codeError * (PDIcode / tau1code)
            oldCodeNco, oldCodeError = codeNco, codeError

            # Save code frequency for current correlation
            trackResults[channelNr]['codeFreq'][loopCnt] = codeFreq
            # Modify code freq based on NCO command
            codeFreq = settings.codeFreqBasis - codeNco

            # --- Record various measures to show in postprocessing ----------------
            tr = trackResults[channelNr]
            tr['dllDiscr'][loopCnt] = codeError
            tr['dllDiscrFilt'][loopCnt] = codeNco
            tr['pllDiscr'][loopCnt] = carrError
            tr['pllDiscrFilt'][loopCnt] = carrNco
            tr['I_E'][loopCnt], tr['I_P'][loopCnt], tr['I_L'][loopCnt] = I_E, I_P, I_L
            tr['Q_E'][loopCnt], tr['Q_P'][loopCnt], tr['Q_L'][loopCnt] = Q_E, Q_P, Q_L

            # --- CNo calculation -------------------------------------------------
            if (loopCnt + 1) % settings.CNo.VSMinterval == 0:
                vsmCnt += 1
                cno_val = cno_vsm(
                    tr['I_P'][loopCnt - settings.CNo.VSMinterval + 1:loopCnt + 1],
                    tr['Q_P'][loopCnt - settings.CNo.VSMinterval + 1:loopCnt + 1],
                    settings.CNo.accTime
                )
                tr['CNo']['VSMValue'][vsmCnt - 1] = cno_val
                tr['CNo']['VSMIndex'][vsmCnt - 1] = loopCnt + 1
                CNo = int(cno_val) if not np.isnan(cno_val) else 'NaN'

        # If we got so far, this means that the tracking was successful
        # Now we only copy status, but it can be updated by a lock detector if implemented
        trackResults[channelNr]['status'] = ch['status']

    return trackResults, channel
