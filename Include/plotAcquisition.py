import matplotlib.pyplot as plt
import numpy as np
import os

def plotAcquisition(acqResults):
    """
    Saves acquisition results plot as .jpg file in the plots directory.
    Args:
        acqResults: Object with .peakMetric (array-like) and .carrFreq (array-like)
    """
    # Create plots directory if it doesn't exist
    plots_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'plots')
    os.makedirs(plots_dir, exist_ok=True)
    
    plt.figure(figsize=(12, 8))
    peakMetric = np.array(acqResults['peakMetric'])
    carrFreq = np.array(acqResults['carrFreq'])
    ax = plt.gca()
    bars = ax.bar(np.arange(1, len(peakMetric)+1), peakMetric, label='Not acquired signals')
    ax.set_title('Acquisition results')
    ax.set_xlabel('PRN number (no bar - SV is not in the acquisition list)')
    ax.set_ylabel('Acquisition Metric')
    ax.set_xlim([0, 33])
    ax.set_xticks(np.arange(1, 33))
    ax.grid(axis='y', which='major')
    acquiredSignals = peakMetric * (carrFreq != 0)
    ax.bar(np.arange(1, len(acquiredSignals)+1), acquiredSignals, color=[0, 0.8, 0], label='Acquired signals')
    ax.legend()
    
    # Save the plot
    filepath = os.path.join(plots_dir, 'acquisition_results.jpg')
    plt.savefig(filepath, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f'   Saved acquisition plot: acquisition_results.jpg')
    print(f'   Acquisition plot saved to: {plots_dir}')
