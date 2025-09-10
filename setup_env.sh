#!/bin/bash

set -e  # Exit immediately on error

cd "$(dirname "$0")"
VENV_DIR="./venv"

# Detect platform and Python
echo "Detected OS: $(uname -s)"
echo "Python version: $(python3 --version)"

# Create venv if missing
if [ ! -d "$VENV_DIR" ]; then
    echo "Creating Python virtual environment..."
    python3 -m venv "$VENV_DIR" || { echo "Failed to create venv"; exit 1; }
    echo "Virtual environment created."
fi

# Activate venv
echo "Activating virtual environment..."
source "$VENV_DIR/bin/activate"

# Upgrade pip tools
echo "Upgrading pip, setuptools, and wheel..."
pip install --upgrade pip setuptools wheel

# Install dependencies from requirements.txt
if [ -f requirements.txt ]; then
    echo "Installing packages from requirements.txt..."
    pip install -r requirements.txt
else
    echo "No requirements.txt found. Skipping package installation."
fi

echo "Setup complete. Virtual environment is ready."
echo "To enter venv, type 'source ./venv/bin/activate'"
echo "Inside of venv, run your code with: python init.py"
