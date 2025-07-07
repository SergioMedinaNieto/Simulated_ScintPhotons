// lar_utils.C
#include "lar_utils.h"
#include "TMath.h"
#include <CLHEP/Vector/ThreeVector.h>
using namespace CLHEP;

// Constants' definitions
double fWion = 23.6 * 1e-6; // MeV
double fWph  = 19.5 * 1e-6; // MeV
double LAr_rho = 1.38434; // g/cm3
double fScintPreScale = 1.52e-1;
double Efield = 0.5;

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

// Quenching parameters
double fQAlpha = 0.72;
double fQProton = 0.868882;

// Parameters' model
bool LArQL = true;
bool fUseBinomialFlucts = true;
bool fUseEllipsModBoxRecomb = true;
bool fUseModBoxRecomb = false;

/* 
Functions' definitions
   Efield : Campo eléctrico (kv/cm)
   dEdx : dEdx (MeV/cm) 
*/

double calcularMedia(std::vector<double> vec) {
    	if (vec.empty()) return 0.0; // Evitar división por cero
    	double suma = std::accumulate(vec.begin(), vec.end(), 0.0);
    	return suma / vec.size();
}

double EfieldAtStep(double x, double y, double z){;
	double shift = +5;
	if (std::abs(x) > 200. + shift || std::abs(y) > 200.+ shift || std::abs(z) > 500.+ shift) return 0.;
  	else 
  		return 0.5; // kV/cm
}

double AngleToEFieldStep(double Efield, double x0, double y0, double z0, double x1, double y1, double z1){
	CLHEP::Hep3Vector v1(x1 - x0, y1 - y0, z1 - z0);
	
	CLHEP::Hep3Vector evec(Efield, 0, 0);
	  if(x0 < 0.)
        	evec.setX(-0.5);
	
	double dotprod = v1.dot(evec);
	double magprod = v1.mag() * evec.mag();
	  if (magprod == 0) return 0.0;
	
	double cosTheta = dotprod / magprod;
	
	if (cosTheta > 1.0) cosTheta = 1.0;
    	if (cosTheta < -1.0) cosTheta = -1.0;
	
	double angle = std::acos(cosTheta);
	
	if (angle > TMath::PiOver2()) { angle = abs(TMath::Pi() - angle); }
	return angle;;	
}

double fBirks(double Efield, double dEdx) {
	return fRecombA / (1. + dEdx * (fRecombk) / Efield);
}

double fModBoxRecomb(double Efield, double dEdx) {
    double Xi = (fModBoxB) * dEdx / Efield;
    return std::log(fModBoxA + Xi) / Xi;
}

double fEscapingEFraction(double dEdx) {
    return fLarqlChi0A / (fLarqlChi0B + std::exp(fLarqlChi0C + fLarqlChi0D * dEdx));
}

double fFieldCorrection(double Efield, double dEdx) {
    return std::exp(-Efield / (fLarqlAlpha * std::log(dEdx) + fLarqlBeta));
}

