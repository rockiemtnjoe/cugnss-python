@echo off
setlocal enabledelayedexpansion

REM
set "SCRIPT_DIR=%~dp0"
cd /d "%SCRIPT_DIR%"
set "VENV_DIR=.\venv"

REM Detect platform and Python
echo Detected OS: Windows
python --version >nul 2>&1
if !errorlevel! neq 0 (
    echo Python is not installed or not in PATH. Please install Python first.
    pause
    exit /b 1
)
echo Python version:
python --version

REM Create venv if missing
if not exist "%VENV_DIR%" (
    echo Creating Python virtual environment...
    python -m venv "%VENV_DIR%"
    if !errorlevel! neq 0 (
        echo Failed to create venv
        pause
        exit /b 1
    )
    echo Virtual environment created.
)

REM Activate venv (Windows)
echo Activating virtual environment...
call "%VENV_DIR%\Scripts\activate.bat"
if !errorlevel! neq 0 (
    echo Failed to activate virtual environment
    pause
    exit /b 1
)

REM Upgrade pip tools
echo Upgrading pip, setuptools, and wheel...
python -m pip install --upgrade pip setuptools wheel
if !errorlevel! neq 0 (
    echo Failed to upgrade pip tools
    pause
    exit /b 1
)

REM Install dependencies from requirements.txt
if exist requirements.txt (
    echo Installing packages from requirements.txt...
    pip install -r requirements.txt
    if !errorlevel! neq 0 (
        echo Failed to install packages from requirements.txt
        pause
        exit /b 1
    )
) else (
    echo No requirements.txt found. Skipping package installation.
)

echo.
echo Setup complete. Virtual environment is ready.
echo To enter venv, type: %VENV_DIR%\Scripts\activate.bat
echo Inside of venv, run your code with: python init.py
echo.
pause
