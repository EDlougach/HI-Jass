@echo off
setlocal
cd /d "%~dp0"

if not exist ".venv\Scripts\python.exe" (
    echo Creating Windows virtual environment...
    py -m venv .venv
    if errorlevel 1 (
        echo Could not create the virtual environment. Install Python 3 and ensure the py launcher is available.
        exit /b 1
    )
)

.venv\Scripts\python.exe -m pip install --upgrade pip
.venv\Scripts\python.exe -m pip install -r requirements.txt
.venv\Scripts\python.exe hi_jass_app.py
