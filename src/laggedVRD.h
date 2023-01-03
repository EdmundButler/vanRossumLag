#pragma once
#include <vector>
typedef std::vector<double> DVEC;

// Entry points from laggedVRD.cpp
bool laggedVRD(const DVEC &sIn, 
           const DVEC &tIn,
           double tau,
           double &P, double &Q, double &Corr,
           double &C);
// Support (assumes tau==1.)               
double VRDcorr(const DVEC &s,
              const DVEC &t,
              double &C);
double VRDnormSQ(const DVEC &x);
double VRDfastCorr(const DVEC &s, const DVEC &t);

// Entry points from WLcorrLag.cpp
struct WeightedTime {double t; double w;};
typedef std::vector<WeightedTime> WTVEC;

//  Main function:
bool WeightedLag(const WTVEC &sIn, 
           const WTVEC &tIn,
           double tau,
           double &P, double &Q, double &Corr,
           double &C);
bool WLcorrLag(const WTVEC &sIn,
               const WTVEC &tIn,
               double tau,
               double &P, double &Q, double &Corr,
               double &C);

double WLcorr(const WTVEC &s,
              const WTVEC &t, 
              double &C);
double WLnormSQ(const WTVEC &x);
double WLfastCorr(const WTVEC &s, const WTVEC &t);

// Test support:
void NoiseData(double Delt, double TotTime, 
               const DVEC& s, DVEC& sMod,
               double p, double q);
void Poisson(double Delt, double TotTime, DVEC &Pvars);

void MakeData(double Delt, double TotTime,
              DVEC &s, DVEC &sMod,
              double p=0., double q=0.);