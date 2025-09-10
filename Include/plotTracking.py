import matplotlib.pyplot as plt
import numpy as np
import os

def plotTracking(channelList, trackResults, settings):
    """
    Saves tracking plots for each channel as .jpg files in the plots directory.
    Args:
        channelList: list of channel indices to plot
        trackResults: list/dict of tracking results per channel
        settings: receiver settings
    """
    channelList = [ch for ch in channelList if 0 <= ch < settings.numberOfChannels]
    
    # Filter to only include channels with tracking data
    valid_channels = [ch for ch in channelList if trackResults[ch]['status'] == 'T']
    
    if not valid_channels:
        print("No valid tracking data to plot")
        return
    
    # Create plots directory if it doesn't exist
    plots_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'plots')
    os.makedirs(plots_dir, exist_ok=True)
    
    # Generate plots for each channel
    for channelNr in valid_channels:
        # Create figure with custom subplot layout
        fig = plt.figure(figsize=(20, 12))
        fig.suptitle(f'Tracking Results for Channel {channelNr} (PRN {trackResults[channelNr]["PRN"]})', fontsize=16)
        
        timeAxisInSeconds = np.arange(1, settings.msToProcess+1) / 1000.0
        
        # Create custom grid layout - 3x3 with some plots spanning 2 columns
        # Top row: scatter plot (1 unit) + nav message (2 units) 
        # Middle row: PLL plots (1 unit each) + correlation results (2 units)
        # Bottom row: DLL plots (1 unit each) + empty space
        
        # Plot 1: Discrete-Time Scatter Plot (top left)
        ax1 = plt.subplot2grid((3, 3), (0, 0), colspan=1, rowspan=1)
        ax1.plot(trackResults[channelNr]['I_P'], trackResults[channelNr]['Q_P'], '.', markersize=2)
        ax1.grid(True)
        ax1.axis('equal')
        ax1.set_title('Discrete-Time Scatter Plot')
        ax1.set_xlabel('I prompt')
        ax1.set_ylabel('Q prompt')
        
        # Plot 2: Navigation Message Bits (top right, spans 2 columns)
        ax2 = plt.subplot2grid((3, 3), (0, 1), colspan=2, rowspan=1)
        ax2.plot(timeAxisInSeconds, trackResults[channelNr]['I_P'], linewidth=1)
        ax2.grid(True)
        ax2.set_title('Bits of the navigation message')
        ax2.set_xlabel('Time (s)')
        ax2.set_ylabel('I prompt')
        
        # Plot 3: Raw PLL Discriminator (middle left)
        ax3 = plt.subplot2grid((3, 3), (1, 0), colspan=1, rowspan=1)
        ax3.plot(timeAxisInSeconds, trackResults[channelNr]['pllDiscr'], 'r', linewidth=1)
        ax3.grid(True)
        ax3.set_xlabel('Time (s)')
        ax3.set_ylabel('Amplitude')
        ax3.set_title('Raw PLL discriminator')
        
        # Plot 4: Correlation Results (middle right, spans 2 columns)
        ax4 = plt.subplot2grid((3, 3), (1, 1), colspan=2, rowspan=1)
        ax4.plot(timeAxisInSeconds, np.sqrt(np.array(trackResults[channelNr]['I_E'])**2 + np.array(trackResults[channelNr]['Q_E'])**2), '-', linewidth=1, label='√(I²E + Q²E)')
        ax4.plot(timeAxisInSeconds, np.sqrt(np.array(trackResults[channelNr]['I_P'])**2 + np.array(trackResults[channelNr]['Q_P'])**2), '-', linewidth=1, label='√(I²P + Q²P)')
        ax4.plot(timeAxisInSeconds, np.sqrt(np.array(trackResults[channelNr]['I_L'])**2 + np.array(trackResults[channelNr]['Q_L'])**2), '-', linewidth=1, label='√(I²L + Q²L)')
        ax4.grid(True)
        ax4.set_title('Correlation results')
        ax4.set_xlabel('Time (s)')
        ax4.set_ylabel('Magnitude')
        ax4.legend(fontsize=8)
        
        # Plot 5: Filtered PLL Discriminator (bottom left)
        ax5 = plt.subplot2grid((3, 3), (2, 0), colspan=1, rowspan=1)
        ax5.plot(timeAxisInSeconds, trackResults[channelNr]['pllDiscrFilt'], 'b', linewidth=1)
        ax5.grid(True)
        ax5.set_xlabel('Time (s)')
        ax5.set_ylabel('Amplitude')
        ax5.set_title('Filtered PLL discriminator')
        
        # Plot 6: Raw DLL Discriminator (bottom middle)
        ax6 = plt.subplot2grid((3, 3), (2, 1), colspan=1, rowspan=1)
        ax6.plot(timeAxisInSeconds, trackResults[channelNr]['dllDiscr'], 'r', linewidth=1)
        ax6.grid(True)
        ax6.set_xlabel('Time (s)')
        ax6.set_ylabel('Amplitude')
        ax6.set_title('Raw DLL discriminator')
        
        # Plot 7: Filtered DLL Discriminator (bottom right)
        ax7 = plt.subplot2grid((3, 3), (2, 2), colspan=1, rowspan=1)
        ax7.plot(timeAxisInSeconds, trackResults[channelNr]['dllDiscrFilt'], 'b', linewidth=1)
        ax7.grid(True)
        ax7.set_xlabel('Time (s)')
        ax7.set_ylabel('Amplitude')
        ax7.set_title('Filtered DLL discriminator')
        
        plt.tight_layout()
        
        # Save the main tracking plot
        filename1 = f'tracking_main_channel_{channelNr}_PRN_{trackResults[channelNr]["PRN"]}.jpg'
        filepath1 = os.path.join(plots_dir, filename1)
        plt.savefig(filepath1, dpi=300, bbox_inches='tight')
        plt.close()
        
        # Separate C/No estimation figure
        fig2, ax_cno = plt.subplots(figsize=(12, 6))
        ax_cno.plot(trackResults[channelNr]['CNo']['VSMValue'], linewidth=1)
        ax_cno.plot(trackResults[channelNr]['CNo']['VSMValue'], 'o', markersize=3)
        ax_cno.set_title(f'CNo Estimation (computed only every 400msec\n(or as specified in initSettings.m)) - Channel {channelNr} PRN {trackResults[channelNr]["PRN"]}')
        ax_cno.set_ylabel('dB-Hz')
        ax_cno.set_xlabel('400msec (or as set in initSettings.m) epoch computation')
        ax_cno.grid(True)
        
        plt.tight_layout()
        
        # Save the C/No plot
        filename2 = f'tracking_cno_channel_{channelNr}_PRN_{trackResults[channelNr]["PRN"]}.jpg'
        filepath2 = os.path.join(plots_dir, filename2)
        plt.savefig(filepath2, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f'   Saved tracking plots: {filename1} and {filename2}')
    
    print(f'   All tracking plots saved to: {plots_dir}')
