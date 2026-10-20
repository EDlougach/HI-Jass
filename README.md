# HI-Jass

Hot Ion Jassby model.

This repository contains a lightweight desktop application for plasma and neutral beam injection (NBI) parameter studies using a compact HotJass-inspired workflow.

The physics lives in its own package, [`hotjass`](https://github.com/EDlougach/hotjass)
(models, User Manual, validation tests), installed by `requirements.txt` at a
pinned version. It is shared with the [HOT-Jass web app](https://github.com/EDlougach/HOT-JASS).

## Features

- Five-tab CustomTkinter interface
- Plasma and NBI parameter input panels
- Live plotting of normalized plasma and NBI profiles
- Simple geometry visualization
- Summary panel for key machine quantities

## Project structure

- `hi_jass_app.py` – desktop GUI entry point
- `requirements.txt` – Python dependencies, including the pinned `hotjass` physics package

To work on the physics together with this app, clone `hotjass` next to this
folder and install it in editable mode: `.venv/bin/pip install -e ../hotjass`.

## Quick start

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python3 hi_jass_app.py
```

### Windows

Install Python 3.10 or newer from python.org, then double-click `run_hi_jass.bat`.
The launcher creates `.venv`, installs the dependencies, and starts the GUI.

From PowerShell or Command Prompt:

```text
cd path\to\HI-Jass
run_hi_jass.bat
```

The GUI must be run in a normal Windows desktop session, not a headless terminal.

See the [User Manual](https://github.com/EDlougach/hotjass/blob/main/User_Manual.md) for the
two-beam slowing-down model, effective fast-ion quantities, density scans, and notation used
by the application.

## GitHub setup

From this folder:

```bash
git init
git add .
git commit -m "Initial HI-Jass app scaffold"
git branch -M main
git remote add origin <your-github-repo-url>
git push -u origin main
```

## Notes

This project is intentionally lightweight and designed as a clean foundation for a more advanced HotJass physics implementation. The current model uses analytic profile shapes to demonstrate the GUI and plotting flow.
