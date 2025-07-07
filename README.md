# LAr Scintillation Photon Simulation Toolkit

![ROOT](https://img.shields.io/badge/ROOT-6.26/10+-FF6D00?logo=root&logoColor=white)
![C++](https://img.shields.io/badge/C++-17-blue?logo=c%2B%2B&logoColor=white)

A simulation tool that models scintillation photon yield in Liquid Argon (LAr) detectors based on energy deposition and electric field effects, with comparison to LArSoft outputs.

## Key Features

- **Three recombination models**:
  - Birks model (default)
  - Modified Box model
  - Ellipsoid Modified Box model
- ROOT-based analysis pipeline
- Direct comparison with LArSoft simulation outputs
- Energy deposition and photon yield histograms

## Repository Structure
LAr_Simulation/
├── code/ # Input .root files
├── images/ # Output plots
├── src/
│ ├── main.C # Main analysis program
│ └── lar_utils.C # Physics utilities (header: lar_utils.h)
└── README.md
