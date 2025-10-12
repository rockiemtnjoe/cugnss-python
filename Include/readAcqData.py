#
#  Read datafile, return a complex numpy vector
#  
import numpy as np
from typing import Tuple
import math
def readAcqData(settings, code_periods = None, skip = None, framing = False) -> np.ndarray:
    """
    read a datafile
    settings is the standard settings object.  We will use:
        fileName
        skipNumberOfBytes
        samplingFreq
        codeFreqBasis
        codeLength
        acqNonCohTime
        acqCoherentInt
    code_periods is the number of 1mS code periods we need
    skip is the number of samples in the datafile to skip
    """
    
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
    if skip == None:
        fid.seek(data_adapt_coeff * settings.skipNumberOfBytes, 0)
    else:
        fid.seek(data_adapt_coeff * skip, 0)
    
    # %% Acquisition ============================================================
    samples_per_code = int(round(settings.codeLength * settings.samplingFreq / settings.codeFreqBasis))
    # At least 42ms of signal are needed for fine frequency estimation

    code_len = (2*settings.acqCoherentInt)*(settings.acqNonCohTime)
    #code_len = max(42, settings.acqNonCohTime + 2)
    
    if code_periods == None:
        num_samples = data_adapt_coeff * code_len * samples_per_code
    else:
        num_samples = code_periods * code_len * samples_per_code
        
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
    fid.close()
    # If the framing flag is set to be true, we fill 1mS of data down each column
    # and there is one column for each code period
    #
    if framing == True:
        data = data.reshape(-1, code_len)
    return (data)