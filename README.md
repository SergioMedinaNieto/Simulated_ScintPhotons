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
// LArQL parameters
double fLarqlChi0A = 0.00338427;
double fLarqlChi0B = -6.57037;
double fLarqlChi0C = 1.88418;
double fLarqlChi0D = 0.000129379;
double fLarqlAlpha = 0.0372;
double fLarqlBeta  = 0.0124;

// Birks parameters
double fRecombA = 0.800;
double fRecombk = 0.0486 / LAr_rho;

// Modified Box Model parameters
double fModBoxA = 0.930;
double fModBoxB = 0.212 / LAr_rho;

// Ellipsoid Modified Box Model parameters
double fEllipsModBoxA =  9.04e-1;
double fEllipsModBoxB =  2.04e-1 / LAr_rho;
double fEllipsModBoxR =  1.25;
```
This file also contains the diferent fuctions based on the model you are using:
```bash
double fBirks(double Efield, double dEdx) {
	return fRecombA / (1. + dEdx * (fRecombk) / Efield);
}

double fModBoxRecomb(double Efield, double dEdx) {
    double Xi = (fModBoxB) * dEdx / Efield;
    return std::log(fModBoxA + Xi) / Xi;
}

double fEllipsModBoxRecomb(double Efield, double dEdx, double phi){
	double B_ellips = fEllipsModBoxB * dEdx / (EfieldStep * std::hypot(std::sin(phi), std::cos(phi) / fEllipsModBoxR));
	return std::log(fEllipsModBoxA + B_ellips) / B_ellips;
}

double fEscapingEFraction(double dEdx) {
    return fLarqlChi0A / (fLarqlChi0B + std::exp(fLarqlChi0C + fLarqlChi0D * dEdx));
}

double fFieldCorrection(double Efield, double dEdx) {
    return std::exp(-Efield / (fLarqlAlpha * std::log(dEdx) + fLarqlBeta));
}
```

## Output file
Example comparison plot:
[Comparison between LArSoft scint photons simulation and our Toy Monte Carlo simulation for a muon (1 GeV) in liquid argon](./images/muon_correlated_NScintPhotonsR.png)



