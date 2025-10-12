import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec 
import numpy as np
from pathlib import Path
import matplotlib.ticker as mtick
from matplotlib.colors import LogNorm, LightSource, PowerNorm


def plotAcqSearch(PRN, settings, results, out_dir=None, dpi=1200, freq_list = None):
    """
    Saves coarse acquisition results plot as .jpg file in the plots directory.
    Args:
        PRN: The PRN of the current satellite; we use it primarily for title & filename
        settings: the settings for the system / sampler
        results: the 2D array of search results; 
                 first index is frequency bins, second index is samples
    """

    nfreq, nsamples = results.shape
    Z = np.asarray(results)

    

    # Shouldn't be complex, but just in case
    if (np.iscomplexobj (Z)):
          Z = np.abs(Z)

    samples_per_code = int(round(settings.samplingFreq / (settings.codeFreqBasis / settings.codeLength)))

    # Create a vector with the frequency steps
    # start_freq = settings.IF - settings.acqSearchBand
    # stop_freq = settings.IF  + settings.acqSearchBand    
    start_freq = settings.acqSearchBand
    stop_freq = -1 * settings.acqSearchBand
    #num_freq = 2 * int(settings.acqSearchBand / settings.acqSearchStep) + 1

    if freq_list == None:
        freq_axis = np.linspace(start_freq, stop_freq, nfreq)
    else:
         assert nfreq == len(freq_list), "Frequency List doesn't match data array size"
         freq_axis = freq_list

 
    samples_per_chip = settings.samplingFreq / settings.codeFreqBasis
    chip_axis = np.arange(nsamples) * 1 / samples_per_chip

    flat_idx = np.argmax(Z)
    peak_freq_idx, peak_chip_idx = np.unravel_index (flat_idx, Z.shape)
 
    # Map indices back to physical units
    peak_freq = freq_axis[peak_freq_idx]
    peak_chip = chip_axis[peak_chip_idx]

    # Gather a subset of information to plot corellation peak
    n_chips = 6    # number of chips we'll plot
    half_width = int (samples_per_chip * n_chips)

    left_idx = max(0, peak_chip_idx - half_width)
    right_idx = min(nsamples, peak_chip_idx + half_width + 1)
    
    x_corr = chip_axis[left_idx:right_idx]
    y_corr = Z[peak_freq_idx, left_idx:right_idx]   # This is what we'll plot

    # Gather the subset for frequency
    n_freq_bins = 7
    left_idx = max(0, peak_freq_idx - n_freq_bins)
    right_idx = min(nfreq, peak_freq_idx + n_freq_bins + 1)
 
    # We'll plot doppler freq in kHz
    x_freq = freq_axis[left_idx:right_idx] * (1 / 1000)
    y_freq = Z[left_idx:right_idx, peak_chip_idx]


    X, Y = np.meshgrid(chip_axis, freq_axis)

    #
    # Now draw the image
    #
    fig = plt.figure(figsize=(8,10))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[2,1])  # A sheet that is 2x2
    ax = fig.add_subplot(gs[0,:], projection="3d")
    #ax.plot3D (X, Y, Z, cmap="coolwarm")

    Zpos = np.clip(Z, np.finfo(float).eps, None)
    vmin = np.percentile(Zpos, 50)
    vmax = np.percentile(Zpos, 98)
    
    norm = PowerNorm(gamma=0.8, vmin=vmin, vmax=vmax)

    #ls = LightSource(azdeg=315, altdeg=45)

    surf = ax.plot_surface (
        X, Y, Zpos, 
        antialiased=False,
        norm=norm,
        cmap="coolwarm",
        linewidth=0,
        rcount=200, ccount=200
        )
    fig.colorbar(surf, ax=ax, shrink=0.6, aspect=12, label="Log-scaled amplitude")
    ax.set_xlabel("Code Phase (chips)")
    ax.set_ylabel("Doppler Frequency (Hz)")
    ax.set_zlabel("Coarse Search - Correlator Output")
    ax.set_title (f"Coarse Acquisition Search - PRN{int(PRN)}")

  
    textstr = (f"Max correlation = {Z[peak_freq_idx, peak_chip_idx]:.2f}\n"
            f"Coarse Doppler = {peak_freq:.0f}, +/- 500 Hz\n"
            f"Chip offset = {peak_chip:.2f} chips \n"
            f"Sample offset = {int(peak_chip_idx)}"
            )

    # Place a text box in upper right corner of the plot
    ax.text2D(0.95, 0.95, textstr,
            transform=ax.transAxes,
            fontsize=10, va="top", ha="right",
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5)
            )
    
    #
    # Draw the peak around maximum corellation
    ax_left = fig.add_subplot(gs[1,0])
    ax_left.plot(x_corr, y_corr, "r", marker = 'o')
    ax_left.set_title ("Time Delay")
    ax_left.set_xlabel ("Chips")

    # Draw the sweep through frequencies
    ax_right = fig.add_subplot(gs[1,1])
    ax_right.plot(x_freq, y_freq, "b", marker = 'o')
    ax_right.set_title ("Doppler")
    ax_right.set_xlabel ("Frequency (kHz)")
    ax_right.xaxis.set_major_formatter(mtick.FormatStrFormatter("%.1f"))

    if out_dir is None:
            out_dir = Path.cwd() / "plots"  # safer than using __file__ in notebooks
    else:
        out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    filepath = out_dir / f"coarse_acquisition_PRN_{int(PRN)}.jpg"
    fig.tight_layout()
    fig.savefig(filepath)
    plt.close(fig)
    return filepath