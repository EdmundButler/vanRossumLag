#include <algorithm>
#include <math.h>
#include "laggedVRD.h"
#include "timing.h"
extern timing sortTime;

// function laggedVRD() - van Rossum distance with lag
//  Inputs: s, t - vectors of increasing spike times
//          tau - time scale
//  Outputs: P, Q - magnitudes of s, t, respectively
//           Corr - correlation of s and t, with optimal lag
//           C - the optimal lag
// returns true if successful
// Model:
//  Let U(t) = Sum(i)[1(t>s(i)) exp((s(i)-t)/tau)]
//      V(t) = Sum(j)[1(t>t(j)) exp((t(i)-t)/tau)]
//  Then
//     P^2 = Integral(U^2(t))
//     Q^2 = Integral(V^2(t))
//     Corr = Integral(U(t) V(t-C))
//     C is the scalar maximizing the expression for Corr
//
//  See vrdLag.pdf for derivation.
bool laggedVRD(const DVEC &sIn, 
           const DVEC &tIn,
           double tau,
           double &P, double &Q, double &Corr,
           double &C)
{
    if (tau <= 0 || sIn.empty() || tIn.empty())
        return false;

    DVEC s(sIn), t(tIn);
    // normalize scale:
    for (double &x : s)
        x /= tau;
    for (double &x : t)
        x /= tau;

    P = sqrt(VRDnormSQ(s));
    Q = sqrt(VRDnormSQ(t));
    Corr = VRDcorr(s,t, C);

    // unscale 
    P *= sqrt(tau);
    Q *= sqrt(tau);
    Corr *= tau;
    C *= tau;

    return true;
}

// function VRDfastCorr() - find dot product of spike vectors
//  Inputs: s, t - the spike vectors
// Returns the correlation without lag (Corr)
//  Implements Algorithm 4 in "vrdLag.pdf"
double VRDfastCorr(const DVEC &s, const DVEC &t)
{
    int M=s.size(), N=t.size();

    // U contains partial j-sums in t[j] > s[i]
    //  Usum accumulates as i decreases
    double U=0., Usum=0.;
    int j=N-1;
    for (int i = M-1; i >= 0; i--)
    {
        if (i<M-1)
            U *= exp(s[i]-s[i+1]);
        for( ; j>=0 && t[j] > s[i]; j--)
            U += exp(s[i]-t[j]);
        Usum += U;
    }

    // V contains partial j-sums in t[j] <= s[i]
    //  Vsum accumulates as i increases
    double V=0., Vsum=0.;
    j = 0;
    for (int i=0; i<M; i++)
    {
        if (i>0)
            V *= exp(s[i-1]-s[i]);
        for( ; j<N && t[j] <= s[i]; j++)
            V += exp(t[j]-s[i]);
        Vsum += V;
    }
    
    return .5*(Usum+Vsum);
}

// function VRDcorr() - find dot product of spike vectors
//  Inputs: s, t - the spike vectors
//  Output: C - the optimal lag           
// Returns the correlation with optimal lag (Corr)
//  Implements Algorithm 2 in "vrdLag.pdf"
double VRDcorr(const DVEC &s,
              const DVEC &t,
              double &C)
{
    int M=s.size(), N=t.size();
    DVEC dels;

    sortTime.startTimer();
    for (int i=0;i<M;i++)
    for (int j=0;j<N;j++)
        dels.push_back(s[i]-t[j]);
    std::sort(dels.begin(), dels.end());
    sortTime.endTimer();
    int MN=dels.size();

    DVEC exps, A({ 1. }), B({ 1. }); // temporary storage
    // need exponentials twice, calculate here
    for (int i=0; i<MN-1; i++)
        exps.push_back(exp(dels[i]-dels[i+1]));
    // Note: B is reversed order compared to document
    for (int i=0; i<MN-1; i++)
    //{ // didnt help
    //    double x = A.back();
    //    double aa = 1. + x - (1.-exps[i])*x;
    //    A.push_back(aa);
    //}
        A.push_back(1. + A.back()*exps[i]);
    for (int i=0; i<MN-1; i++)
        B.push_back(1. + B.back()*exps[MN-i-2]);

    // Need argmax(A[i] + B[MN-i-1])
    double V = 0.;
    int imax=-1;
    for (int i=0; i<MN; i++)
    {
        double Vt = A[i] + B[MN-i-1];
        if (V < Vt)
        {
            imax=i;
            V = Vt;
            C = dels[i]; // output: the optimal lag
        }
    }
    return .5*(V-1.);
}              

// function VRDnormSQ() - find magnitude of spike vector
//  Input: x - increasing times
//  Returns: Integral(U^2(t)), where U(t) = Sum(i)[1(t>x(i)) exp(x(i)-t)]
//   Implements Algorithm 1 in vrdLag.pdf
double VRDnormSQ(const DVEC &x)
{
    int N  = x.size();
    DVEC A({ 1. }); // temporary storage
    
    for (int i=0; i<N-1; i++)
        A.push_back(1. + A.back()*exp(x[i]-x[i+1]));

    double Res = -.5*N;
    for (int i=0; i<N; i++)
        Res += A[i];

    return Res;   
} 
