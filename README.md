# GPS L1 C/A Software Receiver - Python Implementation

A Python-based software-defined radio (SDR) implementation for GPS L1 C/A signal processing. This project is a Python conversion of the CU Multi-GNSS SDR originally developed by  Dennis M. Akos, Darius Plausinaitis, Peter Rinder, and Nicolaj Bertelsen.

## Overview

This software receiver processes raw GPS L1 C/A intermediate frequency (IF) data and performs:
- **Signal Acquisition**: Detects GPS satellites and estimates their carrier frequencies and code phases
- **Signal Tracking**: Tracks carrier phase and code phase for each detected satellite
- **Navigation Data Decoding**: Extracts ephemeris data from navigation messages
- **Position Calculation**: Computes receiver position using pseudorange measurements
- **Visualization**: Generates plots for acquisition, tracking, and navigation results

## Requirements

### System Requirements
- Python 3.7 or higher
- Operating System: Windows, macOS, or Linux

### Python Dependencies
- `numpy` - Numerical computing
- `matplotlib` - Plotting and visualization  
- `scipy` - Scientific computing functions

## Prerequisites

### Python Installation

Before proceeding with the installation, ensure Python is installed on your system and available in your PATH.

#### Check if Python is already installed:
```bash
python --version
# or
python3 --version
```

If Python is installed, you should see output like `Python 3.x.x`. If you get a "command not found" error, you need to install Python.

#### Installing Python:

**For Windows:**
1. Download Python from the [official Python website](https://www.python.org/downloads/windows/)
2. **Important**: During installation, check the box "Add Python to PATH"
3. Follow this detailed tutorial: [How to Install Python on Windows](https://realpython.com/installing-python/#windows)

**For macOS:**
1. **Option 1 - Official Installer**: Download from [python.org](https://www.python.org/downloads/macos/)
2. **Option 2 - Homebrew** (recommended):
   ```bash
   # Install Homebrew if not already installed
   /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
   
   # Install Python
   brew install python
   ```
3. Follow this tutorial: [Installing Python 3 on Mac OS X](https://docs.python-guide.org/starting/install3/osx/)

**For Linux (Ubuntu/Debian):**
```bash
sudo apt update
sudo apt install python3 python3-pip python3-venv
```

#### Verify Installation:
After installation, open a new terminal/command prompt and verify:
```bash
python --version
pip --version
```

Both commands should return version information without errors.

#### Installing Virtual Environment (venv):

The `venv` module is included with Python 3.3+ but may need to be installed separately on some systems:

**For Windows:**
- `venv` is included with Python installations from [python.org](https://www.python.org)
- If missing, reinstall Python and ensure you're using Python 3.7+

**For macOS:**
- `venv` is included with Python installations
- If using Homebrew: `brew install python` includes venv
- If using system Python, you may need: `pip3 install virtualenv`

**For Linux:**
```bash
# Ubuntu/Debian
sudo apt install python3-venv
```

**Verify venv is working:**
```bash
python3 -m venv --help
```
This should display help information without errors.

## Installation

### Option 1: Automated Setup (Recommended)

For macOS/Linux systems:
```bash
# Navigate to project directory
cd GPS_L1CA_Python

# Run setup script
chmod +x setup_env.sh
./setup_env.sh
```

For Windows systems:
```cmd
# Navigate to project directory
cd GPS_L1CA_Python

# Run setup script (double-click or run from command prompt)
setup_env.bat
```

### Option 2: Manual Setup

1. **Clone or download the project**
   ```bash
   git clone <repository-url>
   cd GPS_L1CA_Python
   ```

2. **Create virtual environment**
   ```bash
   python3 -m venv venv
   ```

3. **Activate virtual environment**
   ```bash
   # On macOS/Linux
   source venv/bin/activate
   
   # On Windows
   venv\Scripts\activate
   ```

4. **Install dependencies**
   ```bash
   pip install -r requirements.txt
   ```

## Usage

### Basic Operation

1. **Activate the virtual environment** (if not already active):
   ```bash
   # On macOS/Linux
   source venv/bin/activate
   
   # On Windows (Command Prompt)
   venv\Scripts\activate.bat
   
   # On Windows (PowerShell)
   venv\Scripts\Activate.ps1
   ```

2. **Configure settings** (optional):
   Edit `init_settings.py` to modify:
   - Input data file path (`fileName`)
   - Processing parameters (`msToProcess`, `numberOfChannels`)
   - Acquisition settings (`acqSatelliteList`, `acqSearchBand`)
   - File format settings (`dataType`, `fileType`, `samplingFreq`, `IF`)

3. **Run the software receiver**:
   ```bash
   python init.py
   ```

4. **Follow the interactive prompts**:
   - The software will first probe the input data file
   - Enter `1` to begin GNSS processing or `0` to exit
   - Processing will automatically proceed through acquisition, tracking, and navigation

### Input Data Format

The software expects GPS L1 IF data in binary format. Configure the following in `init_settings.py`:

- **File path**: Set `fileName` to your data file location
- **Data type**: 
  - `'schar'` for 8-bit signed samples
  - `'int16'` for 16-bit signed samples
- **File type**:
  - `1` for real samples (I₀, I₁, I₂, ...)
  - `2` for complex I/Q samples (I₀, Q₀, I₁, Q₁, ...)
- **Sampling frequency**: Set `samplingFreq` in Hz
- **Intermediate frequency**: Set `IF` in Hz

### Output Files

The software generates several output files:

- `acqResults.npy` - Acquisition results (carrier frequencies, code phases, peak metrics)
- `trkResults.npy` - Tracking results (correlator outputs, discriminators, C/N₀)
- `navResults.npy` - Navigation solutions (positions, velocities, timing)
- Various plot files in the `plots/` directory

## Debugging

### VS Code Debugging Setup

To debug the GPS receiver in VS Code, create a `.vscode/launch.json` file in the project root:

```json
{
    "version": "0.2.0",
    "configurations": [
        {
            "name": "Debug GPS Receiver",
            "type": "python",
            "request": "launch",
            "program": "${workspaceFolder}/init.py",
            "console": "integratedTerminal",
            "cwd": "${workspaceFolder}",
            "python": "${workspaceFolder}/venv/bin/python",
            "env": {},
            "args": [],
            "justMyCode": true,
            "stopOnEntry": false
        }
    ]
}
```

**For Windows users**, update the `python` path:
```json
"python": "${workspaceFolder}/venv/Scripts/python.exe"
```

### Selecting the Correct Python Interpreter

#### In VS Code:
1. Open VS Code in the project directory
2. Press `Ctrl+Shift+P` (or `Cmd+Shift+P` on macOS) to open command palette
3. Type "Python: Select Interpreter"
4. Choose the interpreter from your virtual environment:
   - **macOS/Linux**: `./venv/bin/python`
   - **Windows**: `.\venv\Scripts\python.exe`

#### Verify Interpreter Selection:
```bash
# Check which Python is being used
which python  # macOS/Linux
where python  # Windows

# Check Python version and packages
python --version
pip list
```

## Project Structure

```
GPS_L1CA_Python/
├── init.py                 # Main entry point
├── init_settings.py        # Configuration settings
├── requirements.txt        # Python dependencies
├── setup_env.sh/.bat      # Environment setup scripts
├── Common/                 # Shared utility functions
│   ├── calculatePseudoranges.py
│   ├── leastSquarePos.py
│   ├── cart2.py
│   └── ...
├── Include/                # Core processing modules
│   ├── acquisition.py      # Signal acquisition
│   ├── tracking.py         # Signal tracking
│   ├── postNavigation.py   # Navigation processing
│   ├── NAVdecoding.py      # Message decoding
│   └── ...
└── plots/                  # Generated visualization plots
```

## References

This implementation is based on:

1. "A Software-Defined GPS and Galileo Receiver" by Kai Borre et al.
2. CU Multi-GNSS SDR by University of Colorado Boulder
3. GPS Interface Control Documents (ICD-GPS-200)

## License

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

## Authors

- **Original MATLAB Implementation**: Dennis M. Akos, Darius Plausinaitis, Peter Rinder, Nicolaj Bertelsen
- **Python Conversion**: Tyler Schmitz