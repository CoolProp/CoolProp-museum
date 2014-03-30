
#ifndef HELMHOLTZ_H
#define HELMHOLTZ_H

#include <vector>
#include "rapidjson\rapidjson_include.h"

namespace CoolProp{

/// The base class class for the Helmholtz energy terms
/**

Residual Helmholtz Energy Terms:

Term                               | Helmholtz Energy Contribution
----------                         | ------------------------------
ResidualHelmholtzPower             | \f$ \alpha_r=\left\lbrace\begin{array}{cc}\displaystyle\sum_i n_i \delta^{d_i} \tau^{t_i} & l_i=0\\ \displaystyle\sum_i n_i \delta^{d_i} \tau^{t_i} \exp(-\delta^{l_i}) & l_i\neq 0\end{array}\right.\f$
ResidualHelmholtzExponential       | \f$ \alpha_r=\displaystyle\sum_i n_i \delta^{d_i} \tau^{t_i} \exp(-\gamma_i\delta^{l_i}) \f$
ResidualHelmholtzGaussian          | \f$ \alpha_r=\displaystyle\sum_i n_i \delta^{d_i} \tau^{t_i} \exp(-\eta_i(\delta-\epsilon_i)^2-\beta_i(\tau-\gamma_i)^2)\f$
ResidualHelmholtzGERG2008Gaussian  | \f$ \alpha_r=\displaystyle\sum_i n_i \delta^{d_i} \tau^{t_i} \exp(-\eta_i(\delta-\epsilon_i)^2-\beta_i(\delta-\gamma_i))\f$
ResidualHelmholtzLemmon2005        | \f$ \alpha_r=\displaystyle\sum_i n_i \delta^{d_i} \tau^{t_i} \exp(-\delta^{l_i}) \exp(-\tau^{m_i})\f$
ResidualHelmholtzNonAnalytic       | \f$ \begin{array}{c}\alpha_r&=&\displaystyle\sum_i n_i \Delta^{b_i}\delta\psi \\ \Delta & = & \theta^2+B_i[(\delta-1)^2]^{a_i}\\ \theta & = & (1-\tau)+A_i[(\delta-1)^2]^{1/(2\beta_i)}\\ \psi & = & \exp(-C_i(\delta-1)^2-D_i(\tau-1)^2) \end{array}\f$
ResidualHelmholtzSAFTAssociating   | \f$ \alpha_r = am\left(\ln X-\frac{X}{2}+\frac{1}{2}\right); \f$


Ideal-Gas Helmholtz Energy Terms:

Term                                        | Helmholtz Energy Contribution
----------                                  | ------------------------------
IdealHelmholtzLead                          | \f$ \alpha_0 = n_1 + n_2\tau + \ln\delta \f$
IdealHelmholtzEnthalpyEntropyOffset         | \f$ \alpha_0 = \displaystyle\frac{\Delta s}{R_u/M}+\frac{\Delta h}{(R_u/M)T}\tau \f$
IdealHelmholtzLogTau                        | \f$ \alpha_0 = n_1\log\tau \f$
IdealHelmholtzPower                         | \f$ \alpha_0 = \displaystyle\sum_i n_i\tau^{t_i} \f$
IdealHelmholtzPlanckEinstein                | \f$ \alpha_0 = \displaystyle\sum_i n_i\log[1-\exp(-\theta_i\tau)] \f$
IdealHelmholtzPlanckEinstein2               | \f$ \alpha_0 = \displaystyle\sum_i n_i\log[c_i+\exp(\theta_i\tau)] \f$
*/
class BaseHelmholtzTerm{
public:
    BaseHelmholtzTerm(){};
    virtual ~BaseHelmholtzTerm(){};
    /// Returns the base, non-dimensional, Helmholtz energy term (no derivatives) [-]
    /** @param tau Reciprocal reduced temperature where \f$\tau=T_c / T\f$
     *  @param delta Reduced density where \f$\delta = \rho / \rho_c \f$
     */
    virtual long double base(const long double &tau, const long double &delta) throw() = 0;
    /// Returns the first partial derivative of Helmholtz energy term with respect to tau [-]
    /** @param tau Reciprocal reduced temperature where \f$\tau=T_c / T\f$
     *  @param delta Reduced density where \f$\delta = \rho / \rho_c \f$
     */
    virtual long double dTau(const long double &tau, const long double &delta) throw() = 0;
    /// Returns the second partial derivative of Helmholtz energy term with respect to tau [-]
    /** @param tau Reciprocal reduced temperature where \f$\tau=T_c / T\f$
     *  @param delta Reduced density where \f$\delta = \rho / \rho_c \f$
     */ 
    virtual long double dTau2(const long double &tau, const long double &delta) throw() = 0;
    /// Returns the second mixed partial derivative (delta1,dtau1) of Helmholtz energy term with respect to delta and tau [-]
    /** @param tau Reciprocal reduced temperature where \f$\tau=T_c / T\f$
     *  @param delta Reduced density where \f$\delta = \rho / \rho_c \f$
     */
    virtual long double dDelta_dTau(const long double &tau, const long double &delta) throw() = 0;
    /// Returns the first partial derivative of Helmholtz energy term with respect to delta [-]
    /** @param tau Reciprocal reduced temperature where \f$\tau=T_c / T\f$
     *  @param delta Reduced density where \f$\delta = \rho / \rho_c \f$
     */
    virtual long double dDelta(const long double &tau, const long double &delta) throw() = 0;
    /// Returns the second partial derivative of Helmholtz energy term with respect to delta [-]
    /** @param tau Reciprocal reduced temperature where \f$\tau=T_c / T\f$
     *  @param delta Reduced density where \f$\delta = \rho / \rho_c \f$
     */
    virtual long double dDelta2(const long double &tau, const long double &delta) throw() = 0;
    /// Returns the third mixed partial derivative (delta2,dtau1) of Helmholtz energy term with respect to delta and tau [-]
    /** @param tau Reciprocal reduced temperature where \f$\tau=T_c / T\f$
     *  @param delta Reduced density where \f$\delta = \rho / \rho_c \f$
     */
    virtual long double dDelta2_dTau(const long double &tau, const long double &delta) throw() = 0;
    /// Returns the third mixed partial derivative (delta1,dtau2) of Helmholtz energy term with respect to delta and tau [-]
    /** @param tau Reciprocal reduced temperature where \f$\tau=T_c / T\f$
     *  @param delta Reduced density where \f$\delta = \rho / \rho_c \f$
     */
    virtual long double dDelta_dTau2(const long double &tau, const long double &delta) throw() = 0;
    /// Returns the third partial derivative of Helmholtz energy term with respect to tau [-]
    /** @param tau Reciprocal reduced temperature where \f$\tau=T_c / T\f$
     *  @param delta Reduced density where \f$\delta = \rho / \rho_c \f$
     */
    virtual long double dTau3(const long double &tau, const long double &delta) throw() = 0;
    /// Returns the third partial derivative of Helmholtz energy term with respect to delta [-]
    /** @param tau Reciprocal reduced temperature where \f$\tau=T_c / T\f$
     *  @param delta Reduced density where \f$\delta = \rho / \rho_c \f$
     */
    virtual long double dDelta3(const long double &tau, const long double &delta) throw() = 0;
};

// #############################################################################
// #############################################################################
// #############################################################################
//                                RESIDUAL TERMS
// #############################################################################
// #############################################################################
// #############################################################################

struct ResidualHelmholtzPowerElement
{
    long double n,d,t,ld;
    int l;
};
/// Power term
/*!
Term of the form 
\f[ 
\alpha_r=\left\lbrace\begin{array}{cc}\displaystyle\sum_i n_i \delta^{d_i} \tau^{t_i} & l_i=0\\ \displaystyle\sum_i n_i \delta^{d_i} \tau^{t_i} \exp(-\delta^{l_i}) & l_i\neq 0\end{array}\right.
\f]
*/
class ResidualHelmholtzPower : public BaseHelmholtzTerm{
    
public:
    /*std::vector<double> d, ///< The power for the delta terms
                        t, ///< The powers for the tau terms
                        l; ///< The powers for delta in the exp terms*/
    std::vector<long double> s;
    std::size_t N;
    std::vector<ResidualHelmholtzPowerElement> elements;
    // Default Constructor
    ResidualHelmholtzPower(){N = 0;};
    // Constructor
    ResidualHelmholtzPower(const std::vector<long double> &n, const std::vector<long double> &d, 
                           const std::vector<long double> &t, const std::vector<long double> &l)
    {
        N = n.size();
        s.resize(N);
        for (std::size_t i = 0; i < n.size(); ++i)
        {
            ResidualHelmholtzPowerElement el;
            el.n = n[i];
            el.d = d[i];
            el.t = t[i];
            el.ld = l[i];
            el.l = (int)el.ld;
            elements.push_back(el);
        }
    };

    ///< Destructor for the alphar_power class.  No implementation
    ~ResidualHelmholtzPower(){};

    void to_json(rapidjson::Value &el, rapidjson::Document &doc);

    long double base(const long double &tau, const long double &delta);
    long double dDelta(const long double &tau, const long double &delta);
    long double dTau(const long double &tau, const long double &delta);
    long double dDelta2(const long double &tau, const long double &delta);
    long double dDelta_dTau(const long double &tau, const long double &delta);
    long double dTau2(const long double &tau, const long double &delta);
    long double dDelta3(const long double &tau, const long double &delta);
    long double dDelta2_dTau(const long double &tau, const long double &delta);
    long double dDelta_dTau2(const long double &tau, const long double &delta);
    long double dTau3(const long double &tau, const long double &delta);
};

struct ResidualHelmholtzExponentialElement
{
    long double n,d,t,g,l;
};
/**
Term of the form
\f[ \alpha_r=\displaystyle\sum_i n_i \delta^{d_i} \tau^{t_i} \exp(-\gamma_i\delta^{l_i}) \f]
*/
class ResidualHelmholtzExponential : public BaseHelmholtzTerm{

public:
    std::vector<long double> s;
    std::size_t N;
    std::vector<ResidualHelmholtzExponentialElement> elements;
    // Default Constructor
    ResidualHelmholtzExponential(){N = 0;};
    // Constructor
    ResidualHelmholtzExponential(const std::vector<long double> &n, const std::vector<long double> &d, 
                                 const std::vector<long double> &t, const std::vector<long double> &g, 
                                 const std::vector<long double> &l)
    {
        N = n.size();
        s.resize(N);
        for (std::size_t i = 0; i < n.size(); ++i)
        {
            ResidualHelmholtzExponentialElement el;
            el.d = d[i];
            el.t = t[i];
            el.g = g[i];
            el.l = l[i];
            el.n = n[i];
            elements.push_back(el);
        }
    }   

    ///< Destructor for the alphar_power class.  No implementation
    ~ResidualHelmholtzExponential(){};

    void to_json(rapidjson::Value &el, rapidjson::Document &doc);

    long double base(const long double &tau, const long double &delta);
    long double dDelta(const long double &tau, const long double &delta);
    long double dTau(const long double &tau, const long double &delta);
    long double dDelta2(const long double &tau, const long double &delta);
    long double dDelta_dTau(const long double &tau, const long double &delta);
    long double dTau2(const long double &tau, const long double &delta);
    long double dDelta3(const long double &tau, const long double &delta);
    long double dDelta2_dTau(const long double &tau, const long double &delta);
    long double dDelta_dTau2(const long double &tau, const long double &delta);
    long double dTau3(const long double &tau, const long double &delta);
};

struct ResidualHelmholtzGaussianElement
{
    long double n, d, t, eta, epsilon, beta, gamma;
};
class ResidualHelmholtzGaussian : public BaseHelmholtzTerm{

public:
    std::size_t N; ///< The number of terms
    std::vector<ResidualHelmholtzGaussianElement> elements;
    std::vector<long double> s;
    // Default Constructor
    ResidualHelmholtzGaussian(){N = 0;};
    // Constructor
    ResidualHelmholtzGaussian(const std::vector<long double> &n, 
                              const std::vector<long double> &d, 
                              const std::vector<long double> &t, 
                              const std::vector<long double> &eta, 
                              const std::vector<long double> &epsilon,
                              const std::vector<long double> &beta,
                              const std::vector<long double> &gamma
                              )
    { 
        N = n.size(); 
        this->s.resize(N); 
        for (std::size_t i = 0; i < n.size(); ++i)
        {
            ResidualHelmholtzGaussianElement el;
            el.n = n[i];
            el.d = d[i];
            el.t = t[i];
            el.eta = eta[i];
            el.epsilon = epsilon[i];
            el.beta = beta[i];
            el.gamma = gamma[i];
            elements.push_back(el);
        }
    };

    ///< Destructor for the alphar_power class.  No implementation
    ~ResidualHelmholtzGaussian(){};

    void to_json(rapidjson::Value &el, rapidjson::Document &doc);

    long double base(const long double &tau, const long double &delta);
    long double dDelta(const long double &tau, const long double &delta);
    long double dTau(const long double &tau, const long double &delta);
    long double dDelta2(const long double &tau, const long double &delta);
    long double dDelta_dTau(const long double &tau, const long double &delta);
    long double dTau2(const long double &tau, const long double &delta);
    long double dDelta3(const long double &tau, const long double &delta);
    long double dDelta2_dTau(const long double &tau, const long double &delta);
    long double dDelta_dTau2(const long double &tau, const long double &delta);
    long double dTau3(const long double &tau, const long double &delta);
};

class ResidualHelmholtzGERG2008Gaussian : public BaseHelmholtzTerm{

public:
    std::vector<long double> s;
    std::size_t N;
    std::vector<ResidualHelmholtzGaussianElement> elements;
    // Default Constructor
    ResidualHelmholtzGERG2008Gaussian(){N = 0;};
    // Constructor
    ResidualHelmholtzGERG2008Gaussian(const std::vector<long double> &n, 
                                      const std::vector<long double> &d, 
                                      const std::vector<long double> &t, 
                                      const std::vector<long double> &eta, 
                                      const std::vector<long double> &epsilon,
                                      const std::vector<long double> &beta,
                                      const std::vector<long double> &gamma)
    { 
        N = n.size(); 
        s.resize(N); 
        for (std::size_t i = 0; i < n.size(); ++i)
        {
            ResidualHelmholtzGaussianElement el;
            el.n = n[i];
            el.d = d[i];
            el.t = t[i];
            el.eta = eta[i];
            el.epsilon = epsilon[i];
            el.beta = beta[i];
            el.gamma = gamma[i];
            elements.push_back(el);
        }
    };

    ///< Destructor for the alphar_power class.  No implementation
    ~ResidualHelmholtzGERG2008Gaussian(){};

    void to_json(rapidjson::Value &el, rapidjson::Document &doc);

    long double base(const long double &tau, const long double &delta);
    long double dDelta(const long double &tau, const long double &delta);
    long double dTau(const long double &tau, const long double &delta);
    long double dDelta2(const long double &tau, const long double &delta);
    long double dDelta_dTau(const long double &tau, const long double &delta);
    long double dTau2(const long double &tau, const long double &delta);
    long double dDelta3(const long double &tau, const long double &delta){throw NotImplementedError();};
    long double dDelta2_dTau(const long double &tau, const long double &delta){throw NotImplementedError();};
    long double dDelta_dTau2(const long double &tau, const long double &delta){throw NotImplementedError();};
    long double dTau3(const long double &tau, const long double &delta){throw NotImplementedError();};
};

struct ResidualHelmholtzLemmon2005Element{
    long double n, d, t, ld, md;
    int l, m;
};
class ResidualHelmholtzLemmon2005 : public BaseHelmholtzTerm{
public:
    std::size_t N;
    std::vector<long double> s; ///< Summation container
    std::vector<ResidualHelmholtzLemmon2005Element> elements;
    // Default Constructor
    ResidualHelmholtzLemmon2005(){N = 0;};
    // Constructor
    ResidualHelmholtzLemmon2005(const std::vector<long double> &n, 
                                const std::vector<long double> &d, 
                                const std::vector<long double> &t, 
                                const std::vector<long double> &l, 
                                const std::vector<long double> &m)
    {
        N = n.size();
        s.resize(N);
        for (std::size_t i = 0; i < n.size(); ++i)
        {
            ResidualHelmholtzLemmon2005Element el;
            el.n = n[i];
            el.d = d[i];
            el.t = t[i];
            el.ld = l[i];
            el.md = m[i];
            el.l = (int)el.ld;
            el.m = (int)el.md;
            elements.push_back(el);
        }
    };

    ///< Destructor for the alphar_power class.  No implementation
    ~ResidualHelmholtzLemmon2005(){};

    void to_json(rapidjson::Value &el, rapidjson::Document &doc);

    long double base(const long double &tau, const long double &delta);
    long double dDelta(const long double &tau, const long double &delta);
    long double dTau(const long double &tau, const long double &delta);
    long double dDelta2(const long double &tau, const long double &delta);
    long double dDelta_dTau(const long double &tau, const long double &delta);
    long double dTau2(const long double &tau, const long double &delta);
    long double dDelta3(const long double &tau, const long double &delta);
    long double dDelta2_dTau(const long double &tau, const long double &delta);
    long double dDelta_dTau2(const long double &tau, const long double &delta);
    long double dTau3(const long double &tau, const long double &delta);
};

struct ResidualHelmholtzNonAnalyticElement
{
    long double n, a, b, beta, A, B, C, D;
};
class ResidualHelmholtzNonAnalytic : public BaseHelmholtzTerm{

public:
    std::size_t N;
    std::vector<long double> s;
    std::vector<ResidualHelmholtzNonAnalyticElement> elements;
    /// Default Constructor
    ResidualHelmholtzNonAnalytic(){N = 0;};
    /// Destructor. No implementation
    ~ResidualHelmholtzNonAnalytic(){};
    /// Constructor
    ResidualHelmholtzNonAnalytic(const std::vector<long double> &n, 
                                 const std::vector<long double> &a, 
                                 const std::vector<long double> &b, 
                                 const std::vector<long double> &beta, 
                                 const std::vector<long double> &A,
                                 const std::vector<long double> &B,
                                 const std::vector<long double> &C,
                                 const std::vector<long double> &D
                                 )
    {
        N = n.size(); 
        s.resize(N); 
        for (std::size_t i = 0; i < n.size(); ++i)
        {
            ResidualHelmholtzNonAnalyticElement el;
            el.n = n[i];
            el.a = a[i];
            el.b = b[i];
            el.beta = beta[i];
            el.A = A[i];
            el.B = B[i];
            el.C = C[i];
            el.D = D[i];
            elements.push_back(el);
        }
    };

    void to_json(rapidjson::Value &el, rapidjson::Document &doc);

    long double base(const long double &tau, const long double &delta);
    long double dDelta(const long double &tau, const long double &delta);
    long double dTau(const long double &tau, const long double &delta);
    long double dDelta2(const long double &tau, const long double &delta);
    long double dDelta_dTau(const long double &tau, const long double &delta);
    long double dTau2(const long double &tau, const long double &delta);
    long double dDelta3(const long double &tau, const long double &delta);
    long double dDelta2_dTau(const long double &tau, const long double &delta);
    long double dDelta_dTau2(const long double &tau, const long double &delta);
    long double dTau3(const long double &tau, const long double &delta);
};

class ResidualHelmholtzSAFTAssociating : public BaseHelmholtzTerm{
    
protected:
    double a, m,epsilonbar, vbarn, kappabar;

    long double Deltabar(const long double &tau, const long double &delta);
    long double dDeltabar_ddelta__consttau(const long double &tau, const long double &delta);
    long double d2Deltabar_ddelta2__consttau(const long double &tau, const long double &delta);
    long double dDeltabar_dtau__constdelta(const long double &tau, const long double &delta);
    long double d2Deltabar_dtau2__constdelta(const long double &tau, const long double &delta);
    long double d2Deltabar_ddelta_dtau(const long double &tau, const long double &delta);
    long double d3Deltabar_dtau3__constdelta(const long double &tau, const long double &delta);
    long double d3Deltabar_ddelta_dtau2(const long double &tau, const long double &delta);
    long double d3Deltabar_ddelta3__consttau(const long double &tau, const long double &delta);
    long double d3Deltabar_ddelta2_dtau(const long double &tau, const long double &delta);

    long double X(const long double &delta, const long double &Deltabar);
    long double dX_dDeltabar__constdelta(const long double &delta, const long double &Deltabar);
    long double dX_ddelta__constDeltabar(const long double &delta, const long double &Deltabar);
    long double dX_dtau(const long double &tau, const long double &delta);
    long double dX_ddelta(const long double &tau, const long double &delta);
    long double d2X_dtau2(const long double &tau, const long double &delta);
    long double d2X_ddeltadtau(const long double &tau, const long double &delta);
    long double d2X_ddelta2(const long double &tau, const long double &delta);

    long double d3X_dtau3(const long double &tau, const long double &delta);
    long double d3X_ddelta3(const long double &tau, const long double &delta);
    long double d3X_ddeltadtau2(const long double &tau, const long double &delta);
    long double d3X_ddelta2dtau(const long double &tau, const long double &delta);

    long double g(const long double &eta);
    long double dg_deta(const long double &eta);
    long double d2g_deta2(const long double &eta);
    long double d3g_deta3(const long double &eta);
    long double eta(const long double &delta);

public:
    /// Default constructor
    ResidualHelmholtzSAFTAssociating(){ disabled = true; };
    // Constructor
    ResidualHelmholtzSAFTAssociating(double a, double m, double epsilonbar, double vbarn, double kappabar)
        : a(a), m(m), epsilonbar(epsilonbar), vbarn(vbarn), kappabar(kappabar)
    {
        disabled = false;
    };

    bool disabled;

    //Destructor
    ~ResidualHelmholtzSAFTAssociating(){};

    void to_json(rapidjson::Value &el, rapidjson::Document &doc);

    long double base(const long double &tau, const long double &delta);
    long double dDelta(const long double &tau, const long double &delta);
    long double dTau(const long double &tau, const long double &delta);
    long double dDelta2(const long double &tau, const long double &delta);
    long double dDelta_dTau(const long double &tau, const long double &delta);
    long double dTau2(const long double &tau, const long double &delta);
    long double dDelta3(const long double &tau, const long double &delta);
    long double dDelta2_dTau(const long double &tau, const long double &delta);
    long double dDelta_dTau2(const long double &tau, const long double &delta);
    long double dTau3(const long double &tau, const long double &delta);
};

class ResidualHelmholtzContainer
{
    
public:
    ResidualHelmholtzPower Power;
    ResidualHelmholtzExponential Exponential;
    ResidualHelmholtzGaussian Gaussian;
    ResidualHelmholtzLemmon2005 Lemmon2005;
    ResidualHelmholtzNonAnalytic NonAnalytic;
    ResidualHelmholtzSAFTAssociating SAFT;

    long double base(long double tau, long double delta)
    {
        return (Power.base(tau, delta) + Exponential.base(tau, delta)
                +Gaussian.base(tau, delta) + Lemmon2005.base(tau, delta)
                +NonAnalytic.base(tau, delta) + SAFT.base(tau,delta));
    };

    long double dDelta(long double tau, long double delta)
    {
        return (Power.dDelta(tau, delta) + Exponential.dDelta(tau, delta) 
                + Gaussian.dDelta(tau, delta) + Lemmon2005.dDelta(tau, delta) 
                + NonAnalytic.dDelta(tau, delta) + SAFT.dDelta(tau,delta));
    };
    long double dTau(long double tau, long double delta)
    {
        return (Power.dTau(tau, delta) + Exponential.dTau(tau, delta) 
                + Gaussian.dTau(tau, delta) + Lemmon2005.dTau(tau, delta) 
                + NonAnalytic.dTau(tau, delta) + SAFT.dTau(tau,delta));
    };

    long double dDelta2(long double tau, long double delta)
    {
        return (Power.dDelta2(tau, delta) + Exponential.dDelta2(tau, delta)
                +Gaussian.dDelta2(tau, delta) + Lemmon2005.dDelta2(tau, delta)
                +NonAnalytic.dDelta2(tau, delta) + SAFT.dDelta2(tau,delta));
    };
    long double dDelta_dTau(long double tau, long double delta)
    {
        return (Power.dDelta_dTau(tau, delta) + Exponential.dDelta_dTau(tau, delta) 
                + Gaussian.dDelta_dTau(tau, delta) + Lemmon2005.dDelta_dTau(tau, delta) 
                + NonAnalytic.dDelta_dTau(tau, delta) + SAFT.dDelta_dTau(tau,delta));
    };
    long double dTau2(long double tau, long double delta)
    {
        return (Power.dTau2(tau, delta) + Exponential.dTau2(tau, delta) 
                + Gaussian.dTau2(tau, delta) + Lemmon2005.dTau2(tau, delta) 
                + NonAnalytic.dTau2(tau, delta) + SAFT.dTau2(tau,delta));
    };
    
    long double dDelta3(long double tau, long double delta) 
    {
        return (Power.dDelta3(tau, delta) + Exponential.dDelta3(tau, delta)
                +Gaussian.dDelta3(tau, delta) + Lemmon2005.dDelta3(tau, delta)
                +NonAnalytic.dDelta3(tau, delta) + SAFT.dDelta3(tau,delta));
    };
    long double dDelta2_dTau(long double tau, long double delta)
    {
        return (Power.dDelta2_dTau(tau, delta) + Exponential.dDelta2_dTau(tau, delta) 
                + Gaussian.dDelta2_dTau(tau, delta) + Lemmon2005.dDelta2_dTau(tau, delta) 
                + NonAnalytic.dDelta2_dTau(tau, delta) + SAFT.dDelta2_dTau(tau,delta));
    };
    long double dDelta_dTau2(long double tau, long double delta)
    {
        return (Power.dDelta_dTau2(tau, delta) + Exponential.dDelta_dTau2(tau, delta) 
                + Gaussian.dDelta_dTau2(tau, delta) + Lemmon2005.dDelta_dTau2(tau, delta) 
                + NonAnalytic.dDelta_dTau2(tau, delta) + SAFT.dDelta_dTau2(tau,delta));
    };
    long double dTau3(long double tau, long double delta)
    {
        return (Power.dTau3(tau, delta) + Exponential.dTau3(tau, delta)
                +Gaussian.dTau3(tau, delta) + Lemmon2005.dTau3(tau, delta)
                +NonAnalytic.dTau3(tau, delta) + SAFT.dTau3(tau,delta));
    };
};


// #############################################################################
// #############################################################################
// #############################################################################
//                                 IDEAL GAS TERMS
// #############################################################################
// #############################################################################
// #############################################################################

/// The leading term in the EOS used to set the desired reference state
/**
\f[
\alpha_0 = \log(\delta)+a_1+a_2\tau
\f]
*/
class IdealHelmholtzLead : public BaseHelmholtzTerm{

private:
    long double a1, a2;
    bool enabled;
public:
    // Default constructor
    IdealHelmholtzLead(){enabled = false;};

    // Constructor
    IdealHelmholtzLead(const long double a1, const long double a2)
    :a1(a1), a2(a2)
    {enabled = true;};

    //Destructor
    ~IdealHelmholtzLead(){};

    void to_json(rapidjson::Value &el, rapidjson::Document &doc){
        el.AddMember("type","IdealHelmholtzLead",doc.GetAllocator());
        el.AddMember("a1", static_cast<double>(a1), doc.GetAllocator());
        el.AddMember("a2", static_cast<double>(a2), doc.GetAllocator());
    };

    // Term and its derivatives
    long double base(const long double &tau, const long double &delta){
        if (!enabled){return 0.0;}
        return log(delta)+a1+a2*tau;
    };
    long double dDelta(const long double &tau, const long double &delta){
        if (!enabled){return 0.0;}
        return 1.0/delta;
    };
    long double dTau(const long double &tau, const long double &delta){
        if (!enabled){return 0.0;}
        return a2;
    };
    long double dDelta2(const long double &tau, const long double &delta){
        if (!enabled){return 0.0;}
        return -1.0/delta/delta;
    };
    long double dDelta_dTau(const long double &tau, const long double &delta){return 0.0;};
    long double dTau2(const long double &tau, const long double &delta){return 0.0;};
    long double dDelta3(const long double &tau, const long double &delta){
        if (!enabled){return 0.0;}
        return 2/delta/delta/delta;
    };
    long double dDelta2_dTau(const long double &tau, const long double &delta){return 0.0;};
    long double dDelta_dTau2(const long double &tau, const long double &delta){return 0.0;};
    long double dTau3(const long double &tau, const long double &delta){return 0.0;};
};

/// The term in the EOS used to shift the reference state of the fluid
/**
\f[
\alpha_0 = a_1+a_2\tau
\f]
*/
class IdealHelmholtzEnthalpyEntropyOffset : public BaseHelmholtzTerm{
private:
    long double a1,a2; // Use these variables internally
    bool enabled;
public:
    IdealHelmholtzEnthalpyEntropyOffset(){enabled = false;};

    // Constructor
    IdealHelmholtzEnthalpyEntropyOffset(long double a1, long double a2){a1=a1; a2=a2; enabled = true;};

    //Destructor
    ~IdealHelmholtzEnthalpyEntropyOffset(){};

    void to_json(rapidjson::Value &el, rapidjson::Document &doc){
        el.AddMember("type","IdealHelmholtzEnthalpyEntropyOffset",doc.GetAllocator());
        el.AddMember("a1", static_cast<double>(a1), doc.GetAllocator());
        el.AddMember("a2", static_cast<double>(a2), doc.GetAllocator());
    };

    // Term and its derivatives
    long double base(const long double &tau, const long double &delta){
        if (!enabled){return 0.0;}
        return a1+a2*tau;
    };
    long double dDelta(const long double &tau, const long double &delta){return 0.0;};
    long double dTau(const long double &tau, const long double &delta){
        if (!enabled){return 0.0;}
        return a2;
    };
    long double dDelta2(const long double &tau, const long double &delta){return 0.0;};
    long double dDelta_dTau(const long double &tau, const long double &delta){return 0.0;};
    long double dTau2(const long double &tau, const long double &delta){return 0.0;};
    long double dDelta3(const long double &tau, const long double &delta){return 0.0;};
    long double dDelta2_dTau(const long double &tau, const long double &delta){return 0.0;};
    long double dDelta_dTau2(const long double &tau, const long double &delta){return 0.0;};
    long double dTau3(const long double &tau, const long double &delta){return 0.0;};
};


/**
\f[
\alpha_0 = a_1\ln\tau
\f]
*/
class IdealHelmholtzLogTau : public BaseHelmholtzTerm
{
private:
    long double a1;
    bool enabled;
public:

    /// Default constructor
    IdealHelmholtzLogTau(){enabled = false;};

    // Constructor
    IdealHelmholtzLogTau(long double a1){this->a1=a1; enabled = true;};

    //Destructor
    ~IdealHelmholtzLogTau(){};

    void to_json(rapidjson::Value &el, rapidjson::Document &doc){
        el.AddMember("type", "IdealHelmholtzLogTau", doc.GetAllocator());
        el.AddMember("a1", static_cast<double>(a1), doc.GetAllocator());
    };

    // Term and its derivatives
    long double base(const long double &tau, const long double &delta){
        if (!enabled){return 0.0;}
        return a1*log(tau);
    };
    long double dTau(const long double &tau, const long double &delta){
        if (!enabled){return 0.0;}
        return a1/tau;
    };
    long double dTau2(const long double &tau, const long double &delta){
        if (!enabled){return 0.0;}
        return -a1/tau/tau;
    };
    long double dTau3(const long double &tau, const long double &delta){
        if (!enabled){return 0.0;}
        return 2*a1/tau/tau/tau;
    };
    long double dDelta(const long double &tau, const long double &delta){return 0.0;};
    long double dDelta2(const long double &tau, const long double &delta){return 0.0;};
    long double dDelta2_dTau(const long double &tau, const long double &delta){return 0.0;};
    long double dDelta_dTau(const long double &tau, const long double &delta){return 0.0;};
    long double dDelta_dTau2(const long double &tau, const long double &delta){return 0.0;};
    long double dDelta3(const long double &tau, const long double &delta){return 0.0;};
};

/**
\f[
\alpha_0 = \displaystyle\sum_i n_i\tau^{t_i}
\f]
*/
class IdealHelmholtzPower : public BaseHelmholtzTerm{
	
private:
	std::vector<long double> n, t; // Use these variables internally
    std::size_t N;
    bool enabled;
public:
    IdealHelmholtzPower(){enabled = false;};
	// Constructor
	IdealHelmholtzPower(const std::vector<long double> &n, const std::vector<long double> &t)
    :n(n), t(t)
	{
		this->N = n.size();
        enabled = true;
	};

	//Destructor
	~IdealHelmholtzPower(){};

    void to_json(rapidjson::Value &el, rapidjson::Document &doc)
    {
        el.AddMember("type","IdealHelmholtzPower",doc.GetAllocator());
        cpjson::set_long_double_array("n",n,el,doc);
        cpjson::set_long_double_array("t",t,el,doc);
    };

	// Term and its derivatives
	long double base(const long double &tau, const long double &delta){
        if (!enabled){return 0.0;}
		long double s=0; for (std::size_t i = 0; i<N; ++i){s += n[i]*pow(tau, t[i]);} return s;
	};
	long double dTau(const long double &tau, const long double &delta){
        if (!enabled){return 0.0;}
		long double s=0; for (std::size_t i = 0; i<N; ++i){s += n[i]*t[i]*pow(tau, t[i]-1);} return s;
	};
	long double dTau2(const long double &tau, const long double &delta){
        if (!enabled){return 0.0;}
        long double s=0; for (std::size_t i = 0; i<N; ++i){s += n[i]*t[i]*(t[i]-1)*pow(tau, t[i]-2);} return s;
	};
	long double dTau3(const long double &tau, const long double &delta){
        if (!enabled){return 0.0;}
        long double s=0; for (std::size_t i = 0; i<N; ++i){s += n[i]*t[i]*(t[i]-1)*(t[i]-2)*pow(tau, t[i]-3);} return s;
	};
	long double dDelta(const long double &tau, const long double &delta){return 0.0;};
	long double dDelta2(const long double &tau, const long double &delta){return 0.0;};
	long double dDelta2_dTau(const long double &tau, const long double &delta){return 0.0;};
	long double dDelta_dTau(const long double &tau, const long double &delta){return 0.0;};
	long double dDelta_dTau2(const long double &tau, const long double &delta){return 0.0;};
	long double dDelta3(const long double &tau, const long double &delta){return 0.0;};
};

/**
\f[
\alpha_0 = \displaystyle\sum_i n_i\log[1-\exp(-\theta_i\tau)] 
\f]
*/
class IdealHelmholtzPlanckEinstein : public BaseHelmholtzTerm{
	
private:
	std::vector<long double> n,theta; // Use these variables internally
	std::size_t N;
    bool enabled;
public:
    IdealHelmholtzPlanckEinstein(){N = 0; enabled = false;}
	// Constructor with std::vector instances
	IdealHelmholtzPlanckEinstein(std::vector<long double> a, std::vector<long double> theta)
    :n(n), theta(theta)
	{
		N = a.size();
        enabled = false;
	};

	//Destructor
	~IdealHelmholtzPlanckEinstein(){};
  
    void to_json(rapidjson::Value &el, rapidjson::Document &doc)
    {
        el.AddMember("type","IdealHelmholtzPlanckEinstein",doc.GetAllocator());
        cpjson::set_long_double_array("n",n,el,doc);
        cpjson::set_long_double_array("theta",theta,el,doc);
    };

	// Term and its derivatives
	long double base(const long double &tau, const long double &delta){
        if (!enabled){return 0.0;}
        long double s=0; for (std::size_t i=0; i < N; ++i){s += n[i]*log(1.0-exp(-theta[i]*tau));} return s;
    };
    long double dTau(const long double &tau, const long double &delta){
        if (!enabled){return 0.0;}
        long double s=0; for (std::size_t i=0; i < N; ++i){s += n[i]*theta[i]*(1.0/(1.0-exp(-theta[i]*tau))-1.0);} return s;
    };
	long double dTau2(const long double &tau, const long double &delta){
        if (!enabled){return 0.0;}
        long double s=0; for (std::size_t i=0; i < N; ++i){s -= n[i]*pow(theta[i],2)*exp(theta[i]*tau)/pow(1.0-exp(theta[i]*tau),2);} return s;
    };
    long double dTau3(const long double &tau, const long double &delta){
        if (!enabled){return 0.0;}
        long double s=0; for (std::size_t i=0; i < N; ++i){s += n[i]*pow(theta[i],2)*theta[i]*exp(theta[i]*tau)*(exp(theta[i]*tau)+1)/pow(exp(theta[i]*tau)-1,3);} return s;
    };
	long double dDelta(const long double &tau, const long double &delta){return 0.0;};
	long double dDelta2(const long double &tau, const long double &delta){return 0.0;};
	long double dDelta2_dTau(const long double &tau, const long double &delta){return 0.0;};
	long double dDelta_dTau(const long double &tau, const long double &delta){return 0.0;};
	long double dDelta_dTau2(const long double &tau, const long double &delta){return 0.0;};
	long double dDelta3(const long double &tau, const long double &delta){return 0;};
};

/**
\f[ 
\alpha_0 = \displaystyle\sum_i n_i\log[c_i+\exp(\theta_i\tau)] 
\f]
*/
class IdealHelmholtzPlanckEinstein2 : public BaseHelmholtzTerm{
	
private:
	std::vector<long double> n,theta,c; // Use these variables internally
	std::size_t N;
    bool enabled;
public:
    IdealHelmholtzPlanckEinstein2(){N = 0; enabled = false;}
	// Constructor with std::vector instances
	IdealHelmholtzPlanckEinstein2(const std::vector<long double> &n, 
                                  const std::vector<long double> &theta, 
                                  const std::vector<long double> &c)
    :n(n), theta(theta), c(c)
	{
		N = n.size();
        enabled = true;
	};

	//Destructor
	~IdealHelmholtzPlanckEinstein2(){};
  
    void to_json(rapidjson::Value &el, rapidjson::Document &doc)
    {
        el.AddMember("type","IdealHelmholtzPlanckEinstein2",doc.GetAllocator());
        cpjson::set_long_double_array("n",n,el,doc);
        cpjson::set_long_double_array("theta",theta,el,doc);
        cpjson::set_long_double_array("c",c,el,doc);
    };

	// Term and its derivatives
	long double base(const long double &tau, const long double &delta){
        long double s=0; for (std::size_t i=0; i < N; ++i){s += n[i]*log(c[i]+exp(theta[i]*tau));} return s;
    };
    long double dTau(const long double &tau, const long double &delta){
        long double s=0; for (std::size_t i=0; i < N; ++i){s += n[i]*theta[i]*exp(tau*theta[i])/(c[i]+exp(theta[i]*tau));} return s;
    };
	long double dTau2(const long double &tau, const long double &delta){
        long double s=0; for (std::size_t i=0; i < N; ++i){s -= n[i]*pow(theta[i],2)*c[i]*exp(tau*theta[i])/pow(c[i]+exp(tau*theta[i]),2);} return s;
    };
    long double dTau3(const long double &tau, const long double &delta){
        long double s=0; for (std::size_t i=0; i < N; ++i){s += n[i]*pow(theta[i],2)*c[i]*(-theta[i])*exp(theta[i]*tau)*(exp(theta[i]*tau)-c[i])/pow(exp(theta[i]*tau)+c[i],3);} return s;
    };
	long double dDelta(const long double &tau, const long double &delta){return 0.0;};
	long double dDelta2(const long double &tau, const long double &delta){return 0.0;};
	long double dDelta2_dTau(const long double &tau, const long double &delta){return 0.0;};
	long double dDelta_dTau(const long double &tau, const long double &delta){return 0.0;};
	long double dDelta_dTau2(const long double &tau, const long double &delta){return 0.0;};
	long double dDelta3(const long double &tau, const long double &delta){return 0;};
};

class IdealHelmholtzContainer
{
    
public:
    IdealHelmholtzLead Lead;
    IdealHelmholtzEnthalpyEntropyOffset EnthalpyEntropyOffset;
    IdealHelmholtzLogTau LogTau;
    IdealHelmholtzPower Power;
    IdealHelmholtzPlanckEinstein PlanckEinstein;
    IdealHelmholtzPlanckEinstein2 PlanckEinstein2;

    long double base(long double tau, long double delta)
    {
        return (Lead.base(tau, delta) + EnthalpyEntropyOffset.base(tau, delta)
                + LogTau.base(tau, delta) + Power.base(tau, delta) 
                + PlanckEinstein.base(tau, delta) + PlanckEinstein2.base(tau, delta)
                );
    };

    long double dDelta(long double tau, long double delta)
    {
        return 0;
    };
    long double dTau(long double tau, long double delta)
    {
        return 0;
    };
    long double dDelta2(long double tau, long double delta)
    {
        return 0;
    };
    long double dDelta_dTau(long double tau, long double delta)
    {
        return 0;
    };
    long double dTau2(long double tau, long double delta)
    {
        return 0;
    };
    long double dDelta3(long double tau, long double delta) 
    {
        return 0;
    };
    long double dDelta2_dTau(long double tau, long double delta)
    {
        return 0;
    };
    long double dDelta_dTau2(long double tau, long double delta)
    {
        return 0;
    };
    long double dTau3(long double tau, long double delta)
    {
        return 0;
    };
};

};

#endif