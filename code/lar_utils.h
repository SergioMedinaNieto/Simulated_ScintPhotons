// lar_utils.h
#ifndef LAR_UTILS_H
#define LAR_UTILS_H

#include <cmath>
#include <numeric>
#include <vector>
#include <iostream>

extern double calcularMedia(std::vector<double> vec);

extern double EfieldAtStep(double x, double y, double z);

// Constantes
extern double fWion;
extern double fWph;
//extern double Efield;
extern double fScintPreScale;

extern bool LArQL;

extern double Efield;

extern double fLarqlChi0A;
extern double fLarqlChi0B;
extern double fLarqlChi0C;
extern double fLarqlChi0D;

extern double fLarqlAlpha;
extern double fLarqlBeta;

extern double LAr_rho; 

extern bool fUseEllipsModBoxRecomb;
extern double fEllipsModBoxA;
extern double fEllipsModBoxB;
extern double fEllipsModBoxR;

extern bool fUseBinomialFlucts;

extern bool fUseModBoxRecomb;
extern double fRecombA;
extern double fRecombk;

extern double fModBoxA;
extern double fModBoxB;

extern double fQAlpha;
extern double fQProton;

// Funciones
double AngleToEFieldStep(double Efield, double x0, double y0, double z0, double x1, double y1, double z1);
double fModBoxRecomb(double Efield, double dEdx);
double fEscapingEFraction(double dEdx);
double fFieldCorrection(double Efield, double dEdx);
double fBirks(double Efield, double dEdx);

#endif
