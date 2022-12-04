#pragma once
#include <vector>
typedef std::vector<double> DVEC;

// Entry points from laggedVRD.cpp
bool laggedVRD(const DVEC &sIn, 
           const DVEC &tIn,
           double tau,
           double &P, double &Q, double &Corr,
           double &C);
double VRDcorr(const DVEC &s,
              const DVEC &t,
              double &C);
double VRDnormSQ(const DVEC &x);
double VRDfastCorr(const DVEC &s, const DVEC &t);

// Test support:
void NoiseData(double Delt, double TotTime, 
               const DVEC& s, DVEC& sMod,
               double p, double q);

void MakeData(double Delt, double TotTime,
              DVEC &s, DVEC &sMod,
              double p=0., double q=0.);