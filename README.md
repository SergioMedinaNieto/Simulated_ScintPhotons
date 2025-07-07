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

## Requirements

- ROOT (v6.26/10 or later)
- CLHEP libraries
- C++17 compatible compiler

## Installation & Usage

1. **Clone the repository**:
   ```bash
   git clone https://github.com/yourusername/LAr_Simulation.git
   cd LAr_Simulation
2. Prepare input data:
   - Place your .root files in the data/ directory.
   - Example file: muon_correlated.root.
3. Compile and run:
   ```bash
   root -l -q src/main.C
4. Generated plots will be saved in images/

## Configuration
Modify lar_utils.C to select recombination models (Birk model default):
```bash
bool fUseEllipsModBoxRecomb = true;  // Ellipsoid Modified Box
bool fUseModBoxRecomb = false;       // Modified Box
```
You can also see and modify the parameters of each model:
```bash

