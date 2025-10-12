import time
import numpy as np
import cupy as cp
from cupyx.scipy.fft import fft as cp_fft, ifft as cp_ifft
#from cupy.fft import fft as cp_fft, ifft as cp_ifft
from cupyx.scipy.fftpack import get_fft_plan  # plan reuse


def detect_ca_sat_gpu(
    prn: int,
    long_data: np.ndarray,
    settings,
    show_status: bool = False,
    n_coherent: int | None = None,
    n_noncoherent: int | None = None,
) -> CACodeDetectionResult:
    

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

    # ---- PRN template TD/FD
    # Move the PRN onto the GPU, then make a new vector that is
    # copy nCoherent copies and an equal length of zeros
    pn = cp.asarray(make_ca_table(prn, settings), dtype=cp.complex64)
    prn_td = cp.concatenate((cp.tile(pn, n_coherent),
                             cp.zeros(samples_per_code * n_coherent, dtype=cp.complex64)))
    prn_td = cp.ascontiguousarray(prn_td)
    with get_fft_plan(prn_td, value_type='C2C'):
        prn_fd_conj = cp.conj(cp_fft(prn_td))  # (N,)

    # ---- Doppler bins
    freqs = cp.fft.fftfreq(N, t_sample)
    mask = (freqs >= low_f) & (freqs <= high_f)
    F = int(mask.sum().get())   #get() method pulls the sum off the GPU and back to CPU
    initial_shift = (F - 1) // 2
    freq_table_gpu = cp.roll(freqs, initial_shift)[:F].copy()

    # ---- Reshape into B blocks (e.g. B = nNoncoherent) and FFT each (plan once)
    data = x[:total_samples].reshape(B, N)
    data = cp.ascontiguousarray(data)
    with get_fft_plan(data, axes=(1,), value_type='C2C'):
        X_fd = cp_fft(data, axis=1)  # (B, N), contiguous

    # ---- Precompute T_fd for bin 0 and a cheap per-bin update
    # Bin i corresponds to rolling by s = initial_shift - i (mod N).
    # We'll maintain a rolling view of prn_fd_conj using indices (no B×F×N tensors).
    base_idx = cp.arange(N, dtype=cp.int32)
    # start with s0 = initial_shift
    idx = (base_idx - initial_shift) % N
    # preallocate work buffers to avoid re-allocation costs
    conv_fd = cp.empty((B, N), dtype=cp.complex64)
    ifft_out = cp.empty((B, N), dtype=cp.complex64)
    detector_gpu = cp.zeros((F, N), dtype=cp.float32)

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
            power = cp.abs(ifft_out)**2
            detector_gpu[i] = power.sum(axis=0).astype(cp.float32)
            # Update index for next bin: shift right by 1 (equivalent to roll by -1 later)
            idx = (idx + 1) % N

    # Normalize
    detector_gpu *= (1.0 / (B * samples_per_code))

    detector = cp.asnumpy(detector_gpu)
    freq_table = cp.asnumpy(freq_table_gpu)
    return CACodeDetectionResult(True, prn, n_coherent, n_noncoherent, detector, freq_table, time.time()-t0)
