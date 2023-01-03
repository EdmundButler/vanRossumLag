#include <algorithm>
#include <numeric>
#include <math.h>
#include "laggedVRD.h"

#ifdef WANT_TIMING
#include "timing.h"
extern timing sortTime;
#endif

#ifndef SQR
#define SQR(A)((A)*(A))
#endif

// function WeightedLag() - weighted distance with lag
//  Inputs: s, t - vectors of increasing spike times with weights
//          tau - time scale
//  Outputs: P, Q - magnitudes of s, t, respectively
//           Corr - correlation of s and t, with optimal lag
//           C - the optimal lag
// returns true if successful
bool WeightedLag(const WTVEC &sIn, 
           const WTVEC &tIn,
           double tau,
           double &P, double &Q, double &Corr,
           double &C)
{
    if (tau <= 0 || sIn.empty() || tIn.empty())
        return false;

    WTVEC s(sIn), t(tIn);
    // normalize scale:
    for (WeightedTime &x : s)
        x.t /= tau;
    for (WeightedTime &x : t)
        x.t /= tau;

    P = sqrt(WLnormSQ(s));
    Q = sqrt(WLnormSQ(t));
    Corr = WLcorr(s,t, C);

    // unscale 
    P *= sqrt(tau);
    Q *= sqrt(tau);
    Corr *= tau;
    C *= tau;

    return true;
}

// function WLcorrLag() - weighted distance with lag
//  Inputs: s, t - vectors of increasing spike times/weights
//          tau - time scale
//  Outputs: P, Q - magnitudes of s, t, respectively
//           Corr - correlation of s and t, with optimal lag
//           C - the optimal lag
// returns true if successful
// Model:
//  Let U(t) = Sum(i)[1(t>s.t(i)) s.w(i)exp((s.t(i)-t)/tau)]
//      V(t) = Sum(j)[1(t>t.t(j)) t.w(i)exp((t.t(i)-t)/tau)]
//  Then
//     P^2 = Integral(U^2(t))
//     Q^2 = Integral(V^2(t))
//     Corr = Integral(U(t) V(t-C))
//     C is the scalar maximizing the expression for Corr
//
//  See GeneralizedLag.pdf for derivation.
bool WLcorrLag(const WTVEC &sIn, 
           const WTVEC &tIn,
           double tau,
           double &P, double &Q, double &Corr,
           double &C)
{
    if (tau <= 0 || sIn.empty() || tIn.empty())
        return false;

    WTVEC s(sIn), t(tIn);
    // normalize scale:
    for (WeightedTime &x : s)
        x.t /= tau;
    for (WeightedTime &x : t)
        x.t /= tau;

    P = sqrt(WLnormSQ(s));
    Q = sqrt(WLnormSQ(t));
    Corr = WLcorr(s,t, C);

    // unscale 
    P *= sqrt(tau);
    Q *= sqrt(tau);
    Corr *= tau;
    C *= tau;

    return true;
}

// function WLfastCorr() - find dot product of spike vectors
//  Inputs: s, t - the spike vectors
// Returns the correlation without lag (Corr)
//  Implements Algorithm 4 in "GeneralizedLag.pdf"
double WLfastCorr(const WTVEC &s, const WTVEC &t)
{
    int M=s.size(), N=t.size();

    // U contains partial j-sums in t[j] <= s[i]
    //  Vsum accumulates as i increases
    double U=0., Usum=0.;
    int j = 0;
    for (int i=0; i<M; i++)
    {
        if (i>0) // protects against OOB reference to s[]
            U *= exp(s[i-1].t-s[i].t);
        for( ; j<N && t[j].t <= s[i].t; j++)
            U += t[j].w*exp(t[j].t-s[i].t);
        Usum += s[i].w*U;
    }

    // V contains partial j-sums in t[j] > s[i]
    //  Vsum accumulates as i decreases
    double V=0., Vsum=0.;
    j=N-1;
    for (int i = M-1; i >= 0; i--)
    {
        if (i<M-1) // protects against OOB reference to s[]
            V *= exp(s[i].t-s[i+1].t);
        for( ; j>=0 && t[j].t > s[i].t; j--)
            V += t[j].w*exp(s[i].t-t[j].t);
        Vsum += s[i].w*V;
    }
    
    return .5*(Usum+Vsum);
}

// function WLcorr() - find dot product of spike vectors
//  Inputs: s, t - the spike vectors
//  Output: C - the optimal lag           
// Returns the correlation with optimal (ie maximal) lag
//  Implements Algorithm 2 in "GeneralizedLag.pdf"
double WLcorr(const WTVEC &s,
              const WTVEC &t,
              double &C)
{
    /// Sort all differences between s and t in time.
    //   Choose to separate times and weights and sort an index,
    //     to make it easier to determine which sample pair from
    //     s[] and t[] produced the maximum (derive from inds[imax])
    //   Also, perhaps faster than sorting a WTVEC.
    int M=s.size(), N=t.size();
    int MN = M*N;
    DVEC dels, wgts;
#ifdef WANT_TIMING
    sortTime.startTimer();
#endif
    for (int i=0;i<M;i++)
    for (int j=0;j<N;j++)
    {
        dels.push_back(s[i].t-t[j].t);
        wgts.push_back(s[i].w * t[j].w);
    }
    std::vector<int> inds(MN);
    std::iota(inds.begin(), inds.end(), 0); //0,1,...
    class Isrt
    {
        const DVEC &dels;
      public:
        Isrt(const DVEC &delsIn):dels(delsIn){}
        bool operator() (int lhs, int rhs) const 
                       {return dels[lhs]<dels[rhs];}
    }isrt(dels);
    std::sort(inds.begin(), inds.end(), isrt);
#ifdef WANT_TIMING
    sortTime.endTimer();
#endif

    /// Recursions for computing the correlations
    #define D(i) dels[inds[i]]
    #define W(i) wgts[inds[i]]

    DVEC exps, A({ W(0) }), B({ W(MN-1) }); // temporary storage
    // need exponentials twice; calculate here
    for (int i=0; i<MN-1; i++)
        exps.push_back(exp(D(i)-D(i+1)));

    for (int i=0; i<MN-1; i++)
        A.push_back(W(i+1) + A.back()*exps[i]);
    for (int i=MN-2; i>=0; i--)
        B.push_back(W(i) + B.back()*exps[i]);

    /// Find the maximum
    double V = 0.;
    int imax=-1;
    for (int i=0; i<MN; i++)
    {
        double Vt = A[i] + B[MN-i-1];
        if (V < Vt)
        {
            imax=i;
            V = Vt;
            C = D(i); // output: the optimal lag
        }
    }
    return .5*(V-W(imax));
}              

// function WLnormSQ() - find magnitude of spike vector
//  Input: x - increasing times
//  Returns: Integral(U^2(t)), where U(t) = Sum(i)[1(t>x(i)) exp(x(i)-t)]
//   Implements Algorithm 1 in GeneralizedLag.pdf
double WLnormSQ(const WTVEC &x)
{
    int N  = x.size();
    DVEC A({ x[0].w }); // temporary storage
    
    for (int i=0; i<N-1; i++)
        A.push_back(x[i+1].w + A.back()*exp(x[i].t-x[i+1].t));

    double Res = 0;
    for (int i=0; i<N; i++)
        Res += x[i].w*A[i] - .5*SQR(x[i].w);

    return Res;   
} 
