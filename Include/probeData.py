import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import welch

from init_settings import Settings  # Assumes Settings is defined in init_settings.py

def probeData(settings: Settings):
    """
    Plots raw data information: time domain plot, frequency domain plot, and histogram.

    The function can be called as:
        probeData(settings)

    Inputs:
        settings        - receiver settings. Type of data file, sampling
                          frequency and the default filename are specified
                          here.

    --------------------------------------------------------------------------
                                SoftGNSS v3.0

    Copyright (C) Dennis M. Akos
    Written by Darius Plausinaitis and Dennis M. Akos
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
    CVS record:
    $Id: probeData.m,v 1.1.2.7 2006/08/22 13:46:00 dpl Exp $
    _________________________________________________________________________
    """

    fileNameStr = settings.fileName

    # --- Open file and read data ---
    try:
        with open(fileNameStr, 'rb') as f:
            # Move the starting point of processing. Can be used to start the
            # signal processing at any point in the data record (e.g. for long records).
            f.seek(settings.skipNumberOfBytes, 0)
            # Find number of samples per spreading code
            samplesPerCode = int(round(settings.samplingFreq / (settings.codeFreqBasis / settings.codeLength)))
            # Real data: dataAdaptCoeff=1, I/Q data: dataAdaptCoeff=2
            dataAdaptCoeff = 1 if settings.fileType == 1 else 2
            # Read 100ms of signal
            num_samples = dataAdaptCoeff * 100 * samplesPerCode
            # Select data type
            if settings.dataType == 'schar':
                dtype = np.int8
            elif settings.dataType == 'short':
                dtype = np.int16
            elif settings.dataType == 'float':
                dtype = np.float32
            else:
                raise ValueError(f"Unsupported dataType: {settings.dataType}")
            data = np.fromfile(f, dtype=dtype, count=num_samples)
            if data.size < num_samples:
                # The file is too short
                raise ValueError('Could not read enough data from the data file.')
    except Exception as e:
        # Error while opening the data file
        raise RuntimeError(f"Unable to read file {fileNameStr}: {e}")

    # --- Initialization ---------------------------------------------------
    plt.figure(100, figsize=(12, 8))
    plt.clf()
    # Time scale for 5ms
    timeScale = np.arange(0, 5e-3, 1/settings.samplingFreq)
    half_samples = int(round(samplesPerCode/2))

    # --- Time domain plot -------------------------------------------------
    if settings.fileType == 1:
        # Real data
        plt.subplot(2, 2, 3)
        plt.plot(1000 * timeScale[:half_samples], data[:half_samples])
        plt.title('Time domain plot')
        plt.xlabel('Time (ms)')
        plt.ylabel('Amplitude')
        plt.grid(True)
        plt.tight_layout()
    else:
        # I/Q data
        data_cplx = data[::2] + 1j * data[1::2]
        plt.subplot(3, 2, 4)
        plt.plot(1000 * timeScale[:half_samples], np.real(data_cplx[:half_samples]))
        plt.title('Time domain plot (I)')
        plt.xlabel('Time (ms)')
        plt.ylabel('Amplitude')
        plt.grid(True)
        plt.subplot(3, 2, 3)
        plt.plot(1000 * timeScale[:half_samples], np.imag(data_cplx[:half_samples]))
        plt.title('Time domain plot (Q)')
        plt.xlabel('Time (ms)')
        plt.ylabel('Amplitude')
        plt.grid(True)
        plt.tight_layout()

    # --- Frequency domain plot --------------------------------------------
    if settings.fileType == 1:
        # Real Data
        plt.subplot(2, 2, (1, 2))
        f, Pxx = welch(data, fs=settings.samplingFreq, nperseg=32768, noverlap=2048, nfft=32768)
        plt.semilogy(f/1e6, Pxx)
        plt.title('Frequency domain plot')
        plt.xlabel('Frequency (MHz)')
        plt.ylabel('Magnitude')
        plt.grid(True)
    else:
        # I/Q Data
        plt.subplot(3, 2, (1, 2))
        f, Pxx = welch(data_cplx, fs=settings.samplingFreq, nperseg=32768, noverlap=2048, nfft=32768, return_onesided=False)
        # Shift zero freq to center
        f = np.fft.fftshift(f)
        Pxx = np.fft.fftshift(Pxx)
        plt.plot(f/1e6, 10*np.log10(Pxx))
        plt.title('Frequency domain plot')
        plt.xlabel('Frequency (MHz)')
        plt.ylabel('Magnitude (dB)')
        plt.grid(True)

    # --- Histogram --------------------------------------------------------
    if settings.fileType == 1:
        plt.subplot(2, 2, 4)
        plt.hist(data, bins=np.arange(-128, 129))
        dmax = np.max(np.abs(data)) + 1
        plt.axis([-dmax, dmax, None, None])
        plt.title('Histogram')
        plt.xlabel('Bin')
        plt.ylabel('Number in bin')
        plt.grid(True)
    else:
        plt.subplot(3, 2, 6)
        plt.hist(np.real(data_cplx), bins=np.arange(-128, 129))
        dmax = np.max(np.abs(data_cplx)) + 1
        plt.axis([-dmax, dmax, None, None])
        plt.title('Histogram (I)')
        plt.xlabel('Bin')
        plt.ylabel('Number in bin')
        plt.grid(True)
        plt.subplot(3, 2, 5)
        plt.hist(np.imag(data_cplx), bins=np.arange(-128, 129))
        dmax = np.max(np.abs(data_cplx)) + 1
        plt.axis([-dmax, dmax, None, None])
        plt.title('Histogram (Q)')
        plt.xlabel('Bin')
        plt.ylabel('Number in bin')
        plt.grid(True)

    plt.tight_layout()
    plt.show()
