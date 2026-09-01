# HI-Jass

Hot Ion Jassby model.

This repository contains a lightweight desktop application for plasma and neutral beam injection (NBI) parameter studies using a compact HotJass-inspired workflow.

## Features

- Five-tab CustomTkinter interface
- Plasma and NBI parameter input panels
- Live plotting of normalized plasma and NBI profiles
- Simple geometry visualization
- Summary panel for key machine quantities

## Project structure

- `hi_jass_app.py` – desktop GUI entry point
- `hotjass_core.py` – numerical model and parameter definitions
- `requirements.txt` – Python dependencies

## Quick start

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python3 hi_jass_app.py
```

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
