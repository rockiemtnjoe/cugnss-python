# init_settings.py

# -----------------------------------------------------------------------------
#                         CU Multi-GNSS SDR  
# (C) Updated by Yafeng Li, Nagaraj C. Shivaramaiah and Dennis M. Akos
# Based on the original work by Darius Plausinaitis, Peter Rinder, 
# Nicolaj Bertelsen and Dennis M. Akos
# -----------------------------------------------------------------------------
#
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
# -----------------------------------------------------------------------------

from dataclasses import dataclass, field
from typing import List
import numpy as np

@dataclass
class TruePosition:
    # True position of the antenna in UTM system (if known). Otherwise enter
    # all NaN's and mean position will be used as a reference.
    E: float = np.nan
    N: float = np.nan
    U: float = np.nan

@dataclass
class CNoSettings:
    # Accumulation interval in Tracking (in Sec)
    accTime: float = 0.001
    # Accumulation interval for computing VSM C/No (in ms)
    VSMinterval: int = 40

@dataclass
class Settings:
    # =========================================================================
    # Processing settings
    # =========================================================================
    # Number of milliseconds to be processed used 36000 + any transients (see
    # below - in Nav parameters) to ensure nav subframes are provided
    msToProcess: int = 60000  # [ms]

    # Number of channels to be used for signal processing
    numberOfChannels: int = 10

    # Move the starting point of processing. Can be used to start the signal
    # processing at any point in the data record (e.g. for long records). fseek
    # function is used to move the file read point, therefore advance is byte
    # based only.
    skipNumberOfBytes: int = 0

    # =========================================================================
    # Raw signal file name and other parameter
    # =========================================================================
    # This is a "default" name of the data file (signal record) to be used in
    # the post-processing mode
    # fileName: str = '../../../L1_IF20KHz_FS18MHz.bin'
    fileName: str = '../testData/testData.bin'

    # Data type used to store one sample
    dataType: str = 'schar'

    # File Types
    # 1 - 8 bit real samples S0,S1,S2,...
    # 2 - 8 bit I/Q samples I0,Q0,I1,Q1,I2,Q2,...
    fileType: int = 2

    # Intermediate, sampling and code frequencies
    IF: float = 20e3         # [Hz]
    samplingFreq: float = 18e6   # [Hz]
    codeFreqBasis: float = 1.023e6   # [Hz]

    # Define number of chips in a code period
    codeLength: int = 1023

    # =========================================================================
    # Acquisition settings
    # =========================================================================
    # Skips acquisition in the script postProcessing.m if set to 1
    skipAcquisition: int = 0

    # List of satellites to look for. Some satellites can be excluded to speed
    # up acquisition
    acqSatelliteList: List[int] = field(default_factory=lambda: list(range(1, 32)))  # [PRN numbers]

    # Band around IF to search for satellite signal. Depends on max Doppler.
    # It is single sideband, so the whole search band is twice of it.
    acqSearchBand: int = 7000  # [Hz]

    # Non-coherent integration times after 1ms coherent integration
    acqNonCohTime: int = 20    # [ms]

    # Threshold for the signal presence decision rule
    acqThreshold: float = 3.5

    # Frequency search step for coarse acquisition
    acqSearchStep: int = 500   # [Hz]

    # Sampling rate threshold for downsampling 
    resamplingThreshold: float = 8e6  # [Hz]

    # Enable/disable use of downsampling for acquisition
    resamplingflag: int = 0    # 0 - Off; 1 - On

    # =========================================================================
    # Tracking loops settings
    # =========================================================================
    # Code tracking loop parameters
    dllDampingRatio: float = 0.7
    dllNoiseBandwidth: float = 1.5   # [Hz]
    dllCorrelatorSpacing: float = 0.5  # [chips]

    # Carrier tracking loop parameters
    pllDampingRatio: float = 0.7
    pllNoiseBandwidth: float = 20    # [Hz]

    # Integration time for DLL and PLL
    intTime: float = 0.001   # [s]

    # =========================================================================
    # Navigation solution settings
    # =========================================================================
    # Period for calculating pseudoranges and position
    navSolPeriod: int = 500  # [ms]

    # Elevation mask to exclude signals from satellites at low elevation
    elevationMask: int = 5   # [degrees 0 - 90]

    # Enable/disable use of tropospheric correction
    useTropCorr: int = 1     # 0 - Off; 1 - On

    # True position of the antenna in UTM system (if known). Otherwise enter
    # all NaN's and mean position will be used as a reference.
    truePosition: TruePosition = field(default_factory=TruePosition)

    # =========================================================================
    # Plot settings
    # =========================================================================
    # Enable/disable plotting of the tracking results for each channel
    plotTracking: int = 1    # 0 - Off; 1 - On
    plotAcquisition: int = 1
    plotNavigation: int = 1

    # =========================================================================
    # Constants
    # =========================================================================
    c: int = 299792458       # The speed of light, [m/s]
    startOffset: float = 68.802   # [ms] Initial signal travel time

    # =========================================================================
    # CNo Settings
    # =========================================================================
    CNo: CNoSettings = field(default_factory=CNoSettings)

def init_settings() -> Settings:
    """
    Initializes and returns a Settings object with default values.

    Returns:
        settings (Settings): Receiver settings (a structure).
    """
    return Settings()
