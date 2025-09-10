import matplotlib.pyplot as plt
import numpy as np
import os

def plotNavigation(navSolutions, settings):
    """
    Saves navigation plots as .jpg files in the plots directory.
    Args:
        navSolutions: list or dict with navigation solution arrays
        settings: receiver settings (should have .truePosition with E/N/U)
    """
    if navSolutions is None or len(navSolutions) == 0:
        print('plotNavigation: No navigation data to plot.')
        return

    # Handle both list and dict formats - but extract properly for sky plot
    if isinstance(navSolutions, list):
        # If it's a list, assume each entry has the fields we need
        if len(navSolutions) == 0:
            print('plotNavigation: No navigation data to plot.')
            return
        
        # Extract arrays from list of navigation solutions
        E = np.array([sol.get('E', np.nan) for sol in navSolutions if sol is not None])
        N = np.array([sol.get('N', np.nan) for sol in navSolutions if sol is not None])
        U = np.array([sol.get('U', np.nan) for sol in navSolutions if sol is not None])
        longitude = np.array([sol.get('longitude', np.nan) for sol in navSolutions if sol is not None])
        latitude = np.array([sol.get('latitude', np.nan) for sol in navSolutions if sol is not None])
        height = np.array([sol.get('height', np.nan) for sol in navSolutions if sol is not None])
        DOP = np.array([sol.get('DOP', [np.nan, np.nan]) for sol in navSolutions if sol is not None])
        
        # For sky plot, get all unique satellites from the last solution
        if navSolutions and 'az' in navSolutions[-1]:
            last_solution = navSolutions[-1]
            az_data = last_solution['az']
            el_data = last_solution['el']
            prn_data = last_solution['PRN']
            
            # Convert to numpy arrays and handle empty cases
            az = np.array(az_data) if len(az_data) > 0 else np.array([])
            el = np.array(el_data) if len(el_data) > 0 else np.array([])
            PRN = np.array(prn_data) if len(prn_data) > 0 else np.array([])
        else:
            az = np.array([])
            el = np.array([])
            PRN = np.array([])
    else:
        # Original dict format
        E = np.array(navSolutions.get('E', []))
        N = np.array(navSolutions.get('N', []))
        U = np.array(navSolutions.get('U', []))
        longitude = np.array(navSolutions.get('longitude', []))
        latitude = np.array(navSolutions.get('latitude', []))
        height = np.array(navSolutions.get('height', []))
        DOP = navSolutions.get('DOP', np.array([]))
        az = navSolutions.get('az', np.array([]))
        el = navSolutions.get('el', np.array([]))
        PRN = navSolutions.get('PRN', np.array([]))

    if (np.isnan(settings.truePosition.E) or np.isnan(settings.truePosition.N) or np.isnan(settings.truePosition.U)):
        refCoord = {
            'E': np.nanmean(E),
            'N': np.nanmean(N),
            'U': np.nanmean(U)
        }
        refPointLgText = 'Mean Position'
    else:
        refCoord = {
            'E': settings.truePosition.E,
            'N': settings.truePosition.N,
            'U': settings.truePosition.U
        }
        refPointLgText = 'Reference Position'

    # Create plots directory if it doesn't exist
    plots_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'plots')
    os.makedirs(plots_dir, exist_ok=True)
    
    # Create a single figure with multiple subplots
    fig = plt.figure(figsize=(16, 12))
    
    # Plot 1: Coordinate Variations (top)
    ax1 = plt.subplot(3, 2, (1, 2))  # Span across top two columns
    ax1.plot(E - refCoord['E'], label='E', linewidth=1)
    ax1.plot(N - refCoord['N'], label='N', linewidth=1)
    ax1.plot(U - refCoord['U'], label='U', linewidth=1)
    ax1.set_title('Coordinates variations in UTM system')
    ax1.set_xlabel(f'Measurement period: {settings.navSolPeriod}ms')
    ax1.set_ylabel('Variations (m)')
    ax1.grid(True)
    ax1.legend()
    
    # Plot 2: 2D Position Plot (bottom left)
    ax2 = plt.subplot(3, 2, (3, 5))  # Left column, spans middle and bottom
    ax2.plot(E - refCoord['E'], N - refCoord['N'], '+', markersize=4, label='Measurements')
    ax2.scatter(0, 0, c='r', marker='+', s=80, label=refPointLgText)
    ax2.set_title('Positions in UTM system (3D plot)')
    ax2.set_xlabel('East (m)')
    ax2.set_ylabel('North (m)')
    ax2.legend()
    ax2.grid(True)
    ax2.axis('equal')
    
    # Add mean position text
    if len(navSolutions) > 0:
        mean_lat = np.nanmean(latitude)
        mean_lon = np.nanmean(longitude) 
        mean_height = np.nanmean(height)
        
        # Convert latitude to degrees, minutes, seconds
        lat_deg = int(abs(mean_lat))
        lat_min = int((abs(mean_lat) - lat_deg) * 60)
        lat_sec = ((abs(mean_lat) - lat_deg) * 60 - lat_min) * 60
        lat_hemisphere = 'N' if mean_lat >= 0 else 'S'
        
        # Convert longitude to degrees, minutes, seconds
        lon_deg = int(abs(mean_lon))
        lon_min = int((abs(mean_lon) - lon_deg) * 60)
        lon_sec = ((abs(mean_lon) - lon_deg) * 60 - lon_min) * 60
        lon_hemisphere = 'E' if mean_lon >= 0 else 'W'
        
        # Add text box with position info in DMS format
        info_text = f"Measurements\nMean Position\nLat: {lat_deg}°{lat_min:02d}'{lat_sec:05.2f}\"{lat_hemisphere}\nLng: {lon_deg}°{lon_min:02d}'{lon_sec:05.2f}\"{lon_hemisphere}\nHgt: {mean_height:.1f}m"
        ax2.text(0.02, 0.98, info_text, transform=ax2.transAxes, 
                verticalalignment='top', bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    # Plot 3: Sky Plot (bottom right)
    ax3 = plt.subplot(3, 2, (4, 6))  # Right column, spans middle and bottom
    
    # Create a polar-like sky plot with better definition
    angles = np.linspace(0, 2*np.pi, 360)
    
    # Add multiple elevation circles for better definition
    elevation_circles = [15, 30, 45, 60, 75]
    for elev in elevation_circles:
        r = np.cos(np.radians(elev))
        if elev == 30 or elev == 60:
            # Make 30° and 60° circles more prominent
            ax3.plot(r*np.cos(angles), r*np.sin(angles), 'k--', alpha=0.5, linewidth=1)
        else:
            # Lighter circles for other elevations
            ax3.plot(r*np.cos(angles), r*np.sin(angles), 'k--', alpha=0.3, linewidth=0.5)
    
    # Outer circle (horizon) - more prominent
    ax3.plot(np.cos(angles), np.sin(angles), 'k-', linewidth=2)
    
    # Add radial lines for azimuth directions (every 30 degrees)
    for az_deg in range(0, 360, 30):
        az_rad = np.radians(az_deg)
        x_line = [0, np.sin(az_rad)]
        y_line = [0, np.cos(az_rad)]
        if az_deg % 90 == 0:
            # Main compass directions - more prominent
            ax3.plot(x_line, y_line, 'k-', alpha=0.6, linewidth=1)
        else:
            # Secondary directions - lighter
            ax3.plot(x_line, y_line, 'k-', alpha=0.3, linewidth=0.5)
    
    # Add compass directions (azimuth) with better positioning
    ax3.text(0, 1.15, 'N\n0°', ha='center', va='center', fontsize=11, weight='bold')
    ax3.text(1.15, 0, 'E\n90°', ha='center', va='center', fontsize=11, weight='bold')
    ax3.text(0, -1.15, 'S\n180°', ha='center', va='center', fontsize=11, weight='bold')
    ax3.text(-1.15, 0, 'W\n270°', ha='center', va='center', fontsize=11, weight='bold')
    
    # Add intermediate compass directions
    ax3.text(0.82, 0.82, 'NE\n45°', ha='center', va='center', fontsize=9)
    ax3.text(0.82, -0.82, 'SE\n135°', ha='center', va='center', fontsize=9)
    ax3.text(-0.82, -0.82, 'SW\n225°', ha='center', va='center', fontsize=9)
    ax3.text(-0.82, 0.82, 'NW\n315°', ha='center', va='center', fontsize=9)
    
    # Add elevation labels at multiple positions
    ax3.text(0.97, 0, '15°', ha='center', va='center', fontsize=8, alpha=0.7)
    ax3.text(0.87, 0, '30°', ha='center', va='center', fontsize=8, alpha=0.7)
    ax3.text(0.71, 0, '45°', ha='center', va='center', fontsize=8, alpha=0.7)
    ax3.text(0.52, 0, '60°', ha='center', va='center', fontsize=8, alpha=0.7)
    ax3.text(0.26, 0, '75°', ha='center', va='center', fontsize=8, alpha=0.7)
    ax3.text(0.05, 0, '90°', ha='center', va='center', fontsize=8, alpha=0.7)
    
    # Plot satellite positions if available
    if az.size > 0 and el.size > 0:
        # Handle 1D arrays (list of satellites from one measurement)
        for i in range(len(az)):
            if not np.isnan(az[i]) and not np.isnan(el[i]):
                # Convert to sky plot coordinates
                # Elevation angle maps to radius (90° = center, 0° = edge)
                r = np.cos(np.radians(el[i]))
                # Azimuth angle (0° = North, 90° = East)
                x = r * np.sin(np.radians(az[i]))
                y = r * np.cos(np.radians(az[i]))
                
                # Plot satellite position with better visibility
                ax3.plot(x, y, 'ro', markersize=10, markeredgecolor='darkred', markeredgewidth=1)
                
                # Add PRN label with better positioning and visibility
                if len(PRN) > i and not np.isnan(PRN[i]):
                    ax3.text(x, y-0.08, f'{int(PRN[i])}', 
                           ha='center', va='top', fontsize=9, weight='bold', 
                           color='blue', bbox=dict(boxstyle="round,pad=0.2", facecolor="white", alpha=0.8))
    
    # Set limits with a bit more space
    ax3.set_xlim(-1.3, 1.3)
    ax3.set_ylim(-1.3, 1.3)
    ax3.set_aspect('equal')
    ax3.axis('off')
    
    # Add a subtle background
    ax3.add_patch(plt.Circle((0, 0), 1, color='lightgray', alpha=0.1, zorder=-1))
    
    # Add PDOP info
    if DOP.size > 0:
        if DOP.ndim == 2:
            # Get PDOP (typically the second row in DOP matrix)
            pdop_row = min(1, DOP.shape[0] - 1)
            mean_pdop = np.nanmean(DOP[pdop_row, :])
        else:
            mean_pdop = np.nanmean(DOP)
        ax3.set_title(f'Sky plot (mean PDOP: {mean_pdop:.2f})')
    else:
        ax3.set_title('Sky plot')
    
    plt.tight_layout()
    
    # Save the combined navigation plot
    filepath = os.path.join(plots_dir, 'navigation_combined.jpg')
    plt.savefig(filepath, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'   Saved navigation plot: navigation_combined.jpg')
    print(f'   Navigation plot saved to: {plots_dir}')
