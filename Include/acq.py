from __future__ import annotations
from dataclasses import dataclass
from typing import Optional
import time
import numpy as np
from scipy import fft as sfft
import cupy as cp
from cupyx.scipy.fft import fft as cp_fft, ifft as cp_ifft
from cupyx.scipy.fftpack import get_fft_plan  # plan reuse
from Include.makeCaTable import make_ca_table

"""
This is a set of Classes and Functions written by Joe Carey for 
University of Colorado ASEN 6091 in Fall 2025
"""

@dataclass(slots=True)
class CACodeDetectionResult:
    """Container for C/A-code acquisition output."""
    success: bool
    prn: int
    n_coherent: int
    n_noncoherent: int
    detector: Optional[np.ndarray]   # shape: (n_freq_bins, 2*samples_per_code*n_coherent)
    freq_table: Optional[np.ndarray] # shape: (n_freq_bins,)
    runtime: float
    error_message: Optional[str] = None

    # ------------------------------------------------------------
    # Derived information methods
    # ------------------------------------------------------------
    def peak_value(self) -> Optional[float]:
        """
        Return the maximum value in the detector array, or None if unavailable.
        """
        if self.detector is None:
            return None
        return float(np.max(self.detector))

    def peak_location(self) -> Optional[dict]:
        """
        Return the indices and corresponding Doppler frequency of the detector peak.

        Returns
        -------
        dict with keys:
            'freq_index' : int
            'code_index' : int
            'doppler_hz' : float
            'value'      : float
        or None if no detector data is available.
        """
        if self.detector is None or self.freq_table is None:
            return None

        flat_index = int(np.argmax(self.detector))
        freq_index, code_index = np.unravel_index(flat_index, self.detector.shape)
        doppler_hz = float(self.freq_table[freq_index])
        value = float(self.detector[freq_index, code_index])
        return {
            "freq_index": freq_index,
            "code_index": code_index,
            "doppler_hz": doppler_hz,
            "value": value,
        }

    def summary(self) -> str:
        """
        Return a short human-readable summary for quick inspection.
        """
        if not self.success:
            return f"PRN {self.prn}: Acquisition failed ({self.error_message})"
        pk = self.peak_location()
        return (
            f"PRN {self.prn} | "
            f"{self.n_coherent}x{self.n_noncoherent} integration | "
            f"Peak {pk['value']:.3f} at {pk['doppler_hz']:.1f} Hz, "
            f"Index ({pk['freq_index']},{pk['code_index']}) | "
            f"Runtime {self.runtime:.2f}s"
        )

def detect_ca_sat(
    prn: int,
    long_data: np.ndarray,
    settings,
    show_status: bool = False,
    n_coherent: Optional[int] = None,
    n_noncoherent: Optional[int] = None,
) -> CACodeDetectionResult:
    """
    Attempt to detect the GPS C/A code for PRN `prn` in complex baseband `long_data`.

    Returns
    -------
    CACodeDetectionResult
        success, prn, n_coherent, n_noncoherent, detector (2D), freq_table (1D), runtime
    """
    start_time = time.time()
    if show_status:
        print(f"Beginning Detection of PRN {prn}")

    t_sample = 1.0 / settings.samplingFreq
    low_freq = -float(settings.acqSearchBand)
    high_freq = float(settings.acqSearchBand)

    if n_coherent is None:
        n_coherent = int(settings.acqCoherentInt)
    if n_noncoherent is None:
        n_noncoherent = int(settings.acqNonCohTime)

    if show_status:
        print(f"Coherent Integrations: {n_coherent}, Non-Coherent Integrations: {n_noncoherent}")

    # Samples per code and total sample requirement
    samples_per_code = int(round(settings.codeLength * settings.samplingFreq / settings.codeFreqBasis))
    num_codes = 2 * n_coherent * n_noncoherent
    total_samples = samples_per_code * num_codes

    # Ensure enough data
    if len(long_data) < total_samples:
        if show_status:
            print("Insufficient Data to perform acquisition")
            print(f"{total_samples} samples needed, received {len(long_data)}")
        return CACodeDetectionResult(
            success=False,
            prn=prn,
            n_coherent=n_coherent,
            n_noncoherent=n_noncoherent,
            detector=None,
            freq_table=None,
            runtime=time.time() - start_time,
            error_message="Insufficient Data",
        )

    # DC removal & scale
    rl = np.real(long_data)
    im = np.imag(long_data)
    rl -= np.mean(rl)
    im -= np.mean(im)
    max_sample = np.max(np.concatenate((rl, im)))
    long_data = ((0.5 / max_sample) * (rl + 1j * im)).astype(np.complex64)

    # Mix to baseband if needed
    if getattr(settings, "IF", 0.0) != 0.0:
        w = 2.0 * np.pi * t_sample * float(settings.IF)
        phase_points = np.arange(len(long_data), dtype=np.float32) * w
        lo_carrier = np.exp(-1j * phase_points).astype(np.complex64)
        long_data = (long_data * lo_carrier).astype(np.complex64)

    # PRN template (TD and FD)
    pn_code_samples = make_ca_table(prn, settings)
    zero_padding = np.zeros(samples_per_code * n_coherent, dtype=np.complex64)
    prn_template_td = np.concatenate((np.tile(pn_code_samples, n_coherent), zero_padding)).astype(np.complex64)
    prn_template_fd = np.conj(sfft.fft(prn_template_td))

    # Frequency grid and Doppler window
    freqs = np.fft.fftfreq(len(prn_template_td), t_sample)
    freq_mask = (freqs >= low_freq) & (freqs <= high_freq)
    n_freq_bins = int(np.sum(freq_mask))

    # Center near 0 Hz
    initial_shift = (n_freq_bins - 1) // 2
    freq_table = np.roll(freqs, initial_shift)[:n_freq_bins]

    # Reshape input into non-coherent blocks (each: 2 * n_coherent codes)
    n_cols = 2 * samples_per_code * n_coherent
    data = (long_data[:total_samples]).reshape((n_noncoherent, n_cols))

    # Main convolution loop across Doppler bins
    detector = np.zeros((n_freq_bins, n_cols), dtype=np.float32)

    for noncoh_idx in range(n_noncoherent):
        if show_status:
            print(f"Noncoherent Integration {noncoh_idx}")
        coh_detector = np.zeros((n_freq_bins, n_cols), dtype=np.complex64)
        coh_data = data[noncoh_idx, :]

        this_data_fd = sfft.fft(coh_data)
        this_data_fd = np.roll(this_data_fd, initial_shift)

        for i in range(n_freq_bins):
            convolution_fd = this_data_fd * prn_template_fd
            coh_detector[i, :] = sfft.ifft(convolution_fd)
            this_data_fd = np.roll(this_data_fd, -1)

        detector += (np.abs(coh_detector) ** 2).astype(np.float32)

    # Normalize by number of noncoherent sums and samples per code
    detector = (detector / (n_noncoherent * samples_per_code)).astype(np.float32)

    return CACodeDetectionResult(
        success=True,
        prn=prn,
        n_coherent=n_coherent,
        n_noncoherent=n_noncoherent,
        detector=detector,
        freq_table=freq_table,
        runtime=time.time() - start_time,
    )


#
# The following function implements the same basic logic, using the same basic interfaces,
# but accellerated by an NVIDIA GPU, using the CUDA library
#
def detect_ca_sat_gpu(
    prn: int,
    long_data: np.ndarray,
    settings,
    show_status: bool = False,
    n_coherent: int | None = None,
    n_noncoherent: int | None = None,
) -> CACodeDetectionResult:
    #
    #
    t0 = time.time()
    t_sample = 1.0 / settings.samplingFreq
    low_f  = -float(settings.acqSearchBand)
    high_f =  float(settings.acqSearchBand)
    if n_coherent    is None: n_coherent    = int(settings.acqCoherentInt)
    if n_noncoherent is None: n_noncoherent = int(settings.acqNonCohTime)

    samples_per_code = int(round(settings.codeLength * settings.samplingFreq / settings.codeFreqBasis))
    N = 2 * samples_per_code * n_coherent   # samples per coherent vector; IFFT length / columns
    B = n_noncoherent
    total_samples = samples_per_code * (2 * n_coherent * n_noncoherent)

    if len(long_data) < total_samples:
        return CACodeDetectionResult(False, prn, n_coherent, n_noncoherent, None, None, time.time()-t0, "Insufficient Data")

    # ---- Move + preprocess on GPU
    x = cp.asarray(long_data, dtype=cp.complex64)
    rl = x.real; im = x.imag
    rl -= rl.mean(); im -= im.mean()
    scale = cp.maximum(rl.max(), im.max())
    x = ((0.5 / scale) * (rl + 1j*im)).astype(cp.complex64, copy=False)

    IF = float(getattr(settings, "IF", 0.0))
    if IF != 0.0:
        w = cp.float32(2.0 * cp.pi * t_sample * IF)
        n = cp.arange(x.size, dtype=cp.float32)
        x *= cp.exp(-1j * w * n).astype(cp.complex64, copy=False)

    # ---- Compute Doppler bins
    freqs = cp.fft.fftfreq(N, t_sample)
    mask = (freqs >= low_f) & (freqs <= high_f)
    F = int(mask.sum().get())   #get() method pulls the sum off the GPU and back to CPU
    initial_shift = (F - 1) // 2
    freq_table_gpu = cp.roll(freqs, initial_shift)[:F].copy()
    
    # ---- PRN template: extend TD & compute FD
    # Move the PRN onto the GPU, then make a new vector that is
    # copy nCoherent copies and an equal length of zeros
    pn = cp.asarray(make_ca_table(prn, settings), dtype=cp.complex64)
    prn_td = cp.concatenate((cp.tile(pn, n_coherent),
                             cp.zeros(samples_per_code * n_coherent, dtype=cp.complex64)))
    prn_td = cp.ascontiguousarray(prn_td)
    with get_fft_plan(prn_td, value_type='C2C'):
        prn_fd_conj = cp.conj(cp_fft(prn_td))  # (N,)

    # ---- Reshape into B blocks (e.g. B = nNoncoherent) and FFT each (plan once)
    data = x[:total_samples].reshape(B, N)
    data = cp.ascontiguousarray(data)
    with get_fft_plan(data, axes=(1,), value_type='C2C'):
        X_fd = cp_fft(data, axis=1)  # (B, N), contiguous

    # Rather than apply doppler to the sample data, we apply it to the
    # code template.  Since they get multiplied together, this is effectively
    # equivalent.  However, this simplifies things slightly, since the
    # code template remains constant in time (the sample data changes from
    # frame to frame).  We get the doppler shift by changing the way we
    # index into the code template.
    # Bin i corresponds to rolling by s = initial_shift - i (mod N).
    # We'll maintain a rolling view of prn_fd_conj using indices 
    base_idx = cp.arange(N, dtype=cp.int32)
    # start with s0 = initial_shift
    idx = (base_idx + initial_shift) % N
    # preallocate work buffers to avoid re-allocation costs
    conv_fd = cp.empty((B, N), dtype=cp.complex64)
    ifft_out = cp.empty((B, N), dtype=cp.complex64)
    detector_gpu = cp.zeros((F, N), dtype=cp.float32)
    #
    # We normally think about looping over the non-coherent integrations
    # and adding together matrixes several matrixes each of dimension [doppler][sample].
    # It turns out that is not the best approach on a GPU -it is fastest 
    # to loop over the doppler frequencies, and 
    # at the end of the loop perform all the no-coherent integrations 
    # for that doppler frequency.  This is because
    # n_noncoherent < n_frequencies, this minimizes the gather time
    # at the end of the loop.
    #
    # One cuFFT plan for the (B, N) inverse transforms, reused across all bins
    with get_fft_plan(ifft_out, axes=(1,), value_type='C2C'):
        for i in range(F):
            # Gather shifted template spectrum for this Doppler bin
            T_fd = prn_fd_conj[idx]          # (N,)
            # Broadcast multiply (B, N) * (N,) -> (B, N)
            cp.multiply(X_fd, T_fd, out=conv_fd)
            # Batched IFFT over N for all B blocks
            ifft_out[...] = cp_ifft(conv_fd, axis=1)
            # Power & accumulate over B
            # detector row i: sum |ifft_out|^2 across B
            detector_gpu[i] = cp.sum((ifft_out.real * ifft_out.real + ifft_out.imag * ifft_out.imag), axis = 0, dtype = cp.float32)
            #power = cp.abs(ifft_out)**2
            #detector_gpu[i] = power.sum(axis=0).astype(cp.float32)
            # Update index for next bin: shift right by 1 (equivalent to roll by -1 later)
            idx = (idx - 1) % N

    # Normalize
    detector_gpu *= (1.0 / (B * samples_per_code))

    detector = cp.asnumpy(detector_gpu)
    freq_table = cp.asnumpy(freq_table_gpu)
    return CACodeDetectionResult(True, prn, n_coherent, n_noncoherent, detector, freq_table, time.time()-t0)