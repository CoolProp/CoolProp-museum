
#ifndef HELMHOLTZ_H
#define HELMHOLTZ_H

#include <vector>
#include "rapidjson\rapidjson_include.h"

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
IdealHelmholtzLead                          | \f$ \alpha_0 = a_1 + a_2\tau + \ln\delta \f$
IdealHelmholtzEnthalpyEntropyOffset         | \f$ \alpha_0 = \displaystyle\frac{\Delta s}{R_u/M}+\frac{\Delta h}{(R_u/M)T}\tau \f$
IdealHelmholtzLogTau                        | \f$ \alpha_0 = a_1\log\tau \f$
*/
class BaseHelmholtzTerm{
public:
    virtual ~BaseHelmholtzTerm(){};
    /// Returns the base, non-dimensional, Helmholtz energy term (no derivatives) [-]
    /** @param tau Reciprocal reduced temperature where tau=Tc / T
     *  @param delta Reduced density where delta = rho / rhoc 
     */
    virtual double base(double tau,double delta) = 0;
    /// Returns the first partial derivative of Helmholtz energy term with respect to tau [-]
    /** @param tau Reciprocal reduced temperature where tau=Tc / T
     *  @param delta Reduced density where delta = rho / rhoc 
     */
    virtual double dTau(double tau,double delta) = 0;
    /// Returns the second partial derivative of Helmholtz energy term with respect to tau [-]
    /** @param tau Reciprocal reduced temperature where tau=Tc / T
     *  @param delta Reduced density where delta = rho / rhoc 
     */ 
    virtual double dTau2(double tau, double delta) = 0;
    /// Returns the second mixed partial derivative (delta1,dtau1) of Helmholtz energy term with respect to delta and tau [-]
    /** @param tau Reciprocal reduced temperature where tau=Tc / T
     *  @param delta Reduced density where delta = rho / rhoc 
     */
    virtual double dDelta_dTau(double tau, double delta) = 0;
    /// Returns the first partial derivative of Helmholtz energy term with respect to delta [-]
    /** @param tau Reciprocal reduced temperature where tau=Tc / T
     *  @param delta Reduced density where delta = rho / rhoc 
     */
    virtual double dDelta(double tau, double delta) = 0;
    /// Returns the second partial derivative of Helmholtz energy term with respect to delta [-]
    /** @param tau Reciprocal reduced temperature where tau=Tc / T
     *  @param delta Reduced density where delta = rho / rhoc 
     */
    virtual double dDelta2(double tau, double delta) = 0;
    /// Returns the third mixed partial derivative (delta2,dtau1) of Helmholtz energy term with respect to delta and tau [-]
    /** @param tau Reciprocal reduced temperature where tau=Tc / T
     *  @param delta Reduced density where delta = rho / rhoc 
     */
    virtual double dDelta2_dTau(double tau, double delta) = 0;
    /// Returns the third mixed partial derivative (delta1,dtau2) of Helmholtz energy term with respect to delta and tau [-]
    /** @param tau Reciprocal reduced temperature where tau=Tc / T
     *  @param delta Reduced density where delta = rho / rhoc 
     */
    virtual double dDelta_dTau2(double tau, double delta) = 0;
    /// Returns the third partial derivative of Helmholtz energy term with respect to tau [-]
    /** @param tau Reciprocal reduced temperature where tau=Tc / T
     *  @param delta Reduced density where delta = rho / rhoc 
     */
    virtual double dTau3(double tau, double delta) = 0;
    /// Returns the third partial derivative of Helmholtz energy term with respect to delta [-]
    /** @param tau Reciprocal reduced temperature where tau=Tc / T
     *  @param delta Reduced density where delta = rho / rhoc 
     */
    virtual double dDelta3(double tau, double delta) = 0;

    /// Add the data needed for this term into the rapidjson::Value
    /** @param el rapidjson::Value to be filled
     *  @param doc Top-level document that contains the allocator 
     */
    virtual void to_json(rapidjson::Value &el, rapidjson::Document &doc) = 0;
};

class ResidualHelmholtzTerm : public BaseHelmholtzTerm
{
public:
    std::vector<double> n; ///< The coefficients multiplying each term

    enum parameters {iHEA, iHEdA_dDelta, iHEdA_dTau, iHEd2A_dDelta2, iHEd2A_dTau2, iHEd2A_dDelta_dTau, iHEd3A_dDelta3, iHEd3A_dDelta2_dTau, iHEd3A_dDelta_dTau2, iHEd3A_dTau3};

    std::vector<double> dDelta(std::vector<double> tau, std::vector<double> delta);
    std::vector<double> dDelta2(std::vector<double> tau, std::vector<double> delta);
    std::vector<double> dTau(std::vector<double> tau, std::vector<double> delta);
    std::vector<double> dTau2(std::vector<double> tau, std::vector<double> delta);
    std::vector<double> dDelta_dTau(std::vector<double> tau, std::vector<double> delta);

    /// An evaluation function that will evaluate the desired derivative of the A function, where alphar_i = n_i*A_i
    double eval(int key, double tau, double delta);

    double base(double tau, double delta){return eval(iHEA, tau, delta);};
    double dDelta(double tau, double delta){return eval(iHEdA_dDelta, tau, delta);};
    double dTau(double tau, double delta){return eval(iHEdA_dTau, tau, delta);};
    double dDelta2(double tau, double delta){return eval(iHEd2A_dDelta2, tau, delta);};
    double dDelta_dTau(double tau, double delta){return eval(iHEd2A_dDelta_dTau, tau, delta);};
    double dTau2(double tau, double delta){return eval(iHEd2A_dTau2, tau, delta);};
    double dDelta3(double tau, double delta){return eval(iHEd3A_dDelta3, tau, delta);};
    double dDelta2_dTau(double tau, double delta){return eval(iHEd3A_dDelta2_dTau, tau, delta);};
    double dDelta_dTau2(double tau, double delta){return eval(iHEd3A_dDelta_dTau2, tau, delta);};
    double dTau3(double tau, double delta){return eval(iHEd3A_dTau3, tau, delta);};

    // Derivatives for a single term for use in fitter
    /// Returns the base A_i value where alphar = sum(n_i*A_i)
    /** 
     *  @param log_tau Natural logarithm of tau the reciprocal reduced temperature where log_tau = ln(Tc / T)
     *  @param tau  The reciprocal reduced temperature where tau = Tc / T
     *  @param log_delta Natural logarithm of delta the reduced density where delta = rho / rhoc
     *  @param delta Reduced density where delta = rho / rhoc 
     *  @param i The desired index
     */
    virtual double A(double log_tau, double tau, double log_delta, double delta, int i) = 0;
    /// Returns the derivative of the A_i value with respect to delta where alphar = sum(n_i*A_i)
    /** 
     *  @param log_tau Natural logarithm of tau the reciprocal reduced temperature where log_tau = ln(Tc / T)
     *  @param tau  The reciprocal reduced temperature where tau = Tc / T
     *  @param log_delta Natural logarithm of delta the reduced density where delta = rho / rhoc
     *  @param delta Reduced density where delta = rho / rhoc 
     *  @param i The desired index
     */
    virtual double dA_dDelta(double log_tau, double tau, double log_delta, double delta, int i) = 0;
    /// Returns the derivative of the A_i value with respect to tau where alphar = sum(n_i*A_i)
    /** 
     *  @param log_tau Natural logarithm of tau the reciprocal reduced temperature where log_tau = ln(Tc / T)
     *  @param tau  The reciprocal reduced temperature where tau = Tc / T
     *  @param log_delta Natural logarithm of delta the reduced density where delta = rho / rhoc
     *  @param delta Reduced density where delta = rho / rhoc 
     *  @param i The desired index
     */
    virtual double dA_dTau(double log_tau, double tau, double log_delta, double delta, int i) = 0;
    /// Returns the second derivative of the A_i value with respect to tau where alphar = sum(n_i*A_i)
    /** 
     *  @param log_tau Natural logarithm of tau the reciprocal reduced temperature where log_tau = ln(Tc / T)
     *  @param tau  The reciprocal reduced temperature where tau = Tc / T
     *  @param log_delta Natural logarithm of delta the reduced density where delta = rho / rhoc
     *  @param delta Reduced density where delta = rho / rhoc 
     *  @param i The desired index
     */
    virtual double d2A_dTau2(double log_tau, double tau, double log_delta, double delta, int i) = 0;
    /// Returns the second derivative of the A_i value with respect to delta where alphar = sum(n_i*A_i)
    /** 
     *  @param log_tau Natural logarithm of tau the reciprocal reduced temperature where log_tau = ln(Tc / T)
     *  @param tau  The reciprocal reduced temperature where tau = Tc / T
     *  @param log_delta Natural logarithm of delta the reduced density where delta = rho / rhoc
     *  @param delta Reduced density where delta = rho / rhoc 
     *  @param i The desired index
     */
    virtual double d2A_dDelta2(double log_tau, double tau, double log_delta, double delta, int i) = 0;
    /// Returns the second order mixed partial derivative of the A_i value where alphar = sum(n_i*A_i)
    /** 
     *  @param log_tau Natural logarithm of tau the reciprocal reduced temperature where log_tau = ln(Tc / T)
     *  @param tau  The reciprocal reduced temperature where tau = Tc / T
     *  @param log_delta Natural logarithm of delta the reduced density where delta = rho / rhoc
     *  @param delta Reduced density where delta = rho / rhoc 
     *  @param i The desired index
     */
    virtual double d2A_dDelta_dTau(double log_tau, double tau, double log_delta, double delta, int i) = 0;
    /// Returns the third derivative of the A_i value with respect to delta where alphar = sum(n_i*A_i)
    /** 
     *  @param log_tau Natural logarithm of tau the reciprocal reduced temperature where log_tau = ln(Tc / T)
     *  @param tau  The reciprocal reduced temperature where tau = Tc / T
     *  @param log_delta Natural logarithm of delta the reduced density where delta = rho / rhoc
     *  @param delta Reduced density where delta = rho / rhoc 
     *  @param i The desired index
     */
    virtual double d3A_dDelta3(double log_tau, double tau, double log_delta, double delta, int i) = 0;
    /// Returns the third order mixed partial of the A_i value where alphar = sum(n_i*A_i)
    /** 
     *  @param log_tau Natural logarithm of tau the reciprocal reduced temperature where log_tau = ln(Tc / T)
     *  @param tau  The reciprocal reduced temperature where tau = Tc / T
     *  @param log_delta Natural logarithm of delta the reduced density where delta = rho / rhoc
     *  @param delta Reduced density where delta = rho / rhoc 
     *  @param i The desired index
     */
    virtual double d3A_dDelta2_dTau(double log_tau, double tau, double log_delta, double delta, int i) = 0;
    /// Returns the third order mixed partial of the A_i value where alphar = sum(n_i*A_i)
    /** 
     *  @param log_tau Natural logarithm of tau the reciprocal reduced temperature where log_tau = ln(Tc / T)
     *  @param tau  The reciprocal reduced temperature where tau = Tc / T
     *  @param log_delta Natural logarithm of delta the reduced density where delta = rho / rhoc
     *  @param delta Reduced density where delta = rho / rhoc 
     *  @param i The desired index
     */
    virtual double d3A_dDelta_dTau2(double log_tau, double tau, double log_delta, double delta, int i) = 0;
    /// Returns the third derivative of the A_i value with respect to tau where alphar = sum(n_i*A_i)
    /** 
     *  @param log_tau Natural logarithm of tau the reciprocal reduced temperature where log_tau = ln(Tc / T)
     *  @param tau  The reciprocal reduced temperature where tau = Tc / T
     *  @param log_delta Natural logarithm of delta the reduced density where delta = rho / rhoc
     *  @param delta Reduced density where delta = rho / rhoc 
     *  @param i The desired index
     */
    virtual double d3A_dTau3(double log_tau, double tau, double log_delta, double delta, int i) = 0;
};

// #############################################################################
// #############################################################################
// #############################################################################
//                                RESIDUAL TERMS
// #############################################################################
// #############################################################################
// #############################################################################

/// Power term
/*!

Term are of the form 
\f[ \alpha_r=\left\lbrace\begin{array}{cc}\displaystyle\sum_i n_i \delta^{d_i} \tau^{t_i} & l_i=0\\ \displaystyle\sum_i n_i \delta^{d_i} \tau^{t_i} \exp(-\delta^{l_i}) & l_i\neq 0\end{array}\right.\f]

*/
class ResidualHelmholtzPower : public ResidualHelmholtzTerm{
    
public:
    std::vector<double> d, ///< The power for the delta terms
                        t, ///< The powers for the tau terms
                        l; ///< The powers for delta in the exp terms
    // Default Constructor
    ResidualHelmholtzPower(){};
    // Constructor
    ResidualHelmholtzPower(const std::vector<double> &n, const std::vector<double> &d, const std::vector<double> &t, const std::vector<double> &l)
        : d(d), t(t), l(l)
        {this->n = n;};

    ///< Destructor for the alphar_power class.  No implementation
    ~ResidualHelmholtzPower(){};

    void to_json(rapidjson::Value &el, rapidjson::Document &doc);

    /// Derivatives for a single term for use in fitter
    double A(double log_tau, double tau, double log_delta, double delta, int i);
    double dA_dDelta(double log_tau, double tau, double log_delta, double delta, int i);
    double dA_dTau(double log_tau, double tau, double log_delta, double delta, int i);
    double d2A_dTau2(double log_tau, double tau, double log_delta, double delta, int i);
    double d2A_dDelta2(double log_tau, double tau, double log_delta, double delta, int i);
    double d2A_dDelta_dTau(double log_tau, double tau, double log_delta, double delta, int i);
    double d3A_dDelta3(double log_tau, double tau, double log_delta, double delta, int i);
    double d3A_dDelta2_dTau(double log_tau, double tau, double log_delta, double delta, int i);
    double d3A_dDelta_dTau2(double log_tau, double tau, double log_delta, double delta, int i);
    double d3A_dTau3(double log_tau, double tau, double log_delta, double delta, int i);
};

/**
Term of the form
\f[ \alpha_r=\displaystyle\sum_i n_i \delta^{d_i} \tau^{t_i} \exp(-\gamma_i\delta^{l_i}) \f]
*/
class ResidualHelmholtzExponential : public ResidualHelmholtzTerm{

public:
    std::vector<double> d, ///< The power for the delta terms
                        t, ///< The powers for the tau terms
                        g, ///< The coefficient multiplying delta^l_i
                        l; ///< The powers for delta in the exp terms
    // Default Constructor
    ResidualHelmholtzExponential(){};
    // Constructor
    ResidualHelmholtzExponential(const std::vector<double> &n, const std::vector<double> &d, const std::vector<double> &t, const std::vector<double> &g, const std::vector<double> &l)
        : d(d), t(t), g(g), l(l)
        {this->n = n;};

    ///< Destructor for the alphar_power class.  No implementation
    ~ResidualHelmholtzExponential(){};

    void to_json(rapidjson::Value &el, rapidjson::Document &doc);

    /// Derivatives for a single term for use in fitter
    double A(double log_tau, double tau, double log_delta, double delta, int i);
    double dA_dDelta(double log_tau, double tau, double log_delta, double delta, int i);
    double dA_dTau(double log_tau, double tau, double log_delta, double delta, int i);
    double d2A_dTau2(double log_tau, double tau, double log_delta, double delta, int i);
    double d2A_dDelta2(double log_tau, double tau, double log_delta, double delta, int i);
    double d2A_dDelta_dTau(double log_tau, double tau, double log_delta, double delta, int i);
    double d3A_dDelta3(double log_tau, double tau, double log_delta, double delta, int i);
    double d3A_dDelta2_dTau(double log_tau, double tau, double log_delta, double delta, int i);
    double d3A_dDelta_dTau2(double log_tau, double tau, double log_delta, double delta, int i);
    double d3A_dTau3(double log_tau, double tau, double log_delta, double delta, int i);
};


class ResidualHelmholtzGaussian : public ResidualHelmholtzTerm{

public:
    std::vector<double> d, ///< The power for the delta terms
                        t, ///< The powers for the tau terms
                        eta, ///< The coefficient multiplying (delta-epsilon_i)^2
                        epsilon, 
                        beta, ///< The coefficient multiplying (tau-gamma_i)^2
                        gamma;
    // Default Constructor
    ResidualHelmholtzGaussian(){};
    // Constructor
    ResidualHelmholtzGaussian(const std::vector<double> &n, 
                              const std::vector<double> &d, 
                              const std::vector<double> &t, 
                              const std::vector<double> &eta, 
                              const std::vector<double> &epsilon,
                              const std::vector<double> &beta,
                              const std::vector<double> &gamma
                              )
        : d(d), t(t), eta(eta), epsilon(epsilon), beta(beta), gamma(gamma)
        {this->n = n;};

    ///< Destructor for the alphar_power class.  No implementation
    ~ResidualHelmholtzGaussian(){};

    void to_json(rapidjson::Value &el, rapidjson::Document &doc);

    /// Derivatives for a single term for use in fitter
    double A(double log_tau, double tau, double log_delta, double delta, int i);
    double dA_dDelta(double log_tau, double tau, double log_delta, double delta, int i);
    double dA_dTau(double log_tau, double tau, double log_delta, double delta, int i);
    double d2A_dTau2(double log_tau, double tau, double log_delta, double delta, int i);
    double d2A_dDelta2(double log_tau, double tau, double log_delta, double delta, int i);
    double d2A_dDelta_dTau(double log_tau, double tau, double log_delta, double delta, int i);
    double d3A_dDelta3(double log_tau, double tau, double log_delta, double delta, int i);
    double d3A_dDelta2_dTau(double log_tau, double tau, double log_delta, double delta, int i);
    double d3A_dDelta_dTau2(double log_tau, double tau, double log_delta, double delta, int i);
    double d3A_dTau3(double log_tau, double tau, double log_delta, double delta, int i);
};

class ResidualHelmholtzGERG2008Gaussian : public ResidualHelmholtzTerm{

public:
    std::vector<double> d, ///< The power for the delta terms
                        t, ///< The powers for the tau terms
                        eta, ///< The coefficient multiplying (delta-epsilon_i)^2
                        epsilon, 
                        beta, ///< The coefficient multiplying (tau-gamma_i)^2
                        gamma;
    // Default Constructor
    ResidualHelmholtzGERG2008Gaussian(){};
    // Constructor
    ResidualHelmholtzGERG2008Gaussian(const std::vector<double> &n, 
                                      const std::vector<double> &d, 
                                      const std::vector<double> &t, 
                                      const std::vector<double> &eta, 
                                      const std::vector<double> &epsilon,
                                      const std::vector<double> &beta,
                                      const std::vector<double> &gamma
                            )
        : d(d), t(t), eta(eta), epsilon(epsilon), beta(beta), gamma(gamma)
        {this->n = n;};

    ///< Destructor for the alphar_power class.  No implementation
    ~ResidualHelmholtzGERG2008Gaussian(){};

    double psi(double tau, double delta, int i){return exp(-eta[i]*pow(delta-epsilon[i],2)-beta[i]*(delta-gamma[i]));};

    void to_json(rapidjson::Value &el, rapidjson::Document &doc);

    /// Derivatives for a single term for use in fitter
    double A(double log_tau, double tau, double log_delta, double delta, int i);
    double dA_dDelta(double log_tau, double tau, double log_delta, double delta, int i);
    double dA_dTau(double log_tau, double tau, double log_delta, double delta, int i);
    double d2A_dTau2(double log_tau, double tau, double log_delta, double delta, int i);
    double d2A_dDelta2(double log_tau, double tau, double log_delta, double delta, int i);
    double d2A_dDelta_dTau(double log_tau, double tau, double log_delta, double delta, int i);
    double d3A_dDelta3(double log_tau, double tau, double log_delta, double delta, int i){throw ValueError();};;
    double d3A_dDelta2_dTau(double log_tau, double tau, double log_delta, double delta, int i){throw ValueError();};;
    double d3A_dDelta_dTau2(double log_tau, double tau, double log_delta, double delta, int i){throw ValueError();};;
    double d3A_dTau3(double log_tau, double tau, double log_delta, double delta, int i){throw ValueError();};;
};



class ResidualHelmholtzLemmon2005 : public ResidualHelmholtzTerm{
    
public:
    std::vector<double> d, ///< The power for the delta terms
                        t, ///< The powers for the tau terms
                        l, ///< 
                        m; ///<     
    // Default Constructor
    ResidualHelmholtzLemmon2005(){};
    // Constructor
    ResidualHelmholtzLemmon2005(const std::vector<double> &n, 
                                const std::vector<double> &d, 
                                const std::vector<double> &t, 
                                const std::vector<double> &l, 
                                const std::vector<double> &m
                                )
        : d(d), t(t), l(l), m(m)
        {this->n = n;};

    ///< Destructor for the alphar_power class.  No implementation
    ~ResidualHelmholtzLemmon2005(){};

    void to_json(rapidjson::Value &el, rapidjson::Document &doc);

    /// Derivatives for a single term for use in fitter
    double A(double log_tau, double tau, double log_delta, double delta, int i);
    double dA_dDelta(double log_tau, double tau, double log_delta, double delta, int i);
    double dA_dTau(double log_tau, double tau, double log_delta, double delta, int i);
    double d2A_dTau2(double log_tau, double tau, double log_delta, double delta, int i);
    double d2A_dDelta2(double log_tau, double tau, double log_delta, double delta, int i);
    double d2A_dDelta_dTau(double log_tau, double tau, double log_delta, double delta, int i);
    double d3A_dDelta3(double log_tau, double tau, double log_delta, double delta, int i);
    double d3A_dDelta2_dTau(double log_tau, double tau, double log_delta, double delta, int i);
    double d3A_dDelta_dTau2(double log_tau, double tau, double log_delta, double delta, int i);
    double d3A_dTau3(double log_tau, double tau, double log_delta, double delta, int i);
};

class ResidualHelmholtzNonAnalytic : public ResidualHelmholtzTerm{

public:
    std::vector<double> a, ///< 
                        b, ///< 
                        beta,
                        _A,
                        B,
                        C,
                        D;
    // Default Constructor
    ResidualHelmholtzNonAnalytic(){};
    // Constructor
    ResidualHelmholtzNonAnalytic(const std::vector<double> &n, 
                                 const std::vector<double> &a, 
                                 const std::vector<double> &b, 
                                 const std::vector<double> &beta, 
                                 const std::vector<double> &A,
                                 const std::vector<double> &B,
                                 const std::vector<double> &C,
                                 const std::vector<double> &D
                                 )
        : a(a), b(b), beta(beta), _A(A), B(B), C(C), D(D)
        {this->n = n;};

    ///< Destructor for the alphar_power class.  No implementation
    ~ResidualHelmholtzNonAnalytic(){};

    void to_json(rapidjson::Value &el, rapidjson::Document &doc);

    /// Derivatives for a single term for use in fitter
    double A(double log_tau, double tau, double log_delta, double delta, int i);
    double dA_dDelta(double log_tau, double tau, double log_delta, double delta, int i);
    double dA_dTau(double log_tau, double tau, double log_delta, double delta, int i);
    double d2A_dTau2(double log_tau, double tau, double log_delta, double delta, int i);
    double d2A_dDelta2(double log_tau, double tau, double log_delta, double delta, int i);
    double d2A_dDelta_dTau(double log_tau, double tau, double log_delta, double delta, int i);
    double d3A_dDelta3(double log_tau, double tau, double log_delta, double delta, int i);
    double d3A_dDelta2_dTau(double log_tau, double tau, double log_delta, double delta, int i);
    double d3A_dDelta_dTau2(double log_tau, double tau, double log_delta, double delta, int i);
    double d3A_dTau3(double log_tau, double tau, double log_delta, double delta, int i);
};

class ResidualHelmholtzSAFTAssociating : public ResidualHelmholtzTerm{
    
protected:
    double a, m,epsilonbar, vbarn, kappabar;

    double Deltabar(double tau, double delta);
    double dDeltabar_ddelta__consttau(double tau, double delta);
    double d2Deltabar_ddelta2__consttau(double tau, double delta);
    double dDeltabar_dtau__constdelta(double tau, double delta);
    double d2Deltabar_dtau2__constdelta(double tau, double delta);
    double d2Deltabar_ddelta_dtau(double tau, double delta);
    double d3Deltabar_dtau3__constdelta(double tau, double delta);
    double d3Deltabar_ddelta_dtau2(double tau, double delta);
    double d3Deltabar_ddelta3__consttau(double tau, double delta);
    double d3Deltabar_ddelta2_dtau(double tau, double delta);

    double X(double delta, double Deltabar);
    double dX_dDeltabar__constdelta(double delta, double Deltabar);
    double dX_ddelta__constDeltabar(double delta, double Deltabar);
    double dX_dtau(double tau, double delta);
    double dX_ddelta(double tau, double delta);
    double d2X_dtau2(double tau, double delta);
    double d2X_ddeltadtau(double tau, double delta);
    double d2X_ddelta2(double tau, double delta);

    double d3X_dtau3(double tau, double delta);
    double d3X_ddelta3(double tau, double delta);
    double d3X_ddeltadtau2(double tau, double delta);
    double d3X_ddelta2dtau(double tau, double delta);

    double g(double eta);
    double dg_deta(double eta);
    double d2g_deta2(double eta);
    double d3g_deta3(double eta);
    double eta(double delta);

public:
    // Constructor
    ResidualHelmholtzSAFTAssociating(double a, double m, double epsilonbar, double vbarn, double kappabar)
        : a(a), m(m), epsilonbar(epsilonbar), vbarn(vbarn), kappabar(kappabar)
        {n = std::vector<double>(1,1);};

    //Destructor
    ~ResidualHelmholtzSAFTAssociating(){};

    void to_json(rapidjson::Value &el, rapidjson::Document &doc);

    double A(double log_tau, double tau, double log_delta, double delta, int i);
    double dA_dDelta(double log_tau, double tau, double log_delta, double delta, int i);
    double dA_dTau(double log_tau, double tau, double log_delta, double delta, int i);
    
    double d2A_dDelta2(double log_tau, double tau, double log_delta, double delta, int i);
    double d2A_dDelta_dTau(double log_tau, double tau, double log_delta, double delta, int i);
    double d2A_dTau2(double log_tau, double tau, double log_delta, double delta, int i);
    
    double d3A_dDelta3(double log_tau, double tau, double log_delta, double delta, int i);
    double d3A_dDelta2_dTau(double log_tau, double tau, double log_delta, double delta, int i);
    double d3A_dDelta_dTau2(double log_tau, double tau, double log_delta, double delta, int i);
    double d3A_dTau3(double log_tau, double tau, double log_delta, double delta, int i);
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
    double a1, a2;
public:
    // Constructor
    IdealHelmholtzLead(double a1, double a2){a1=a1; a2=a2;};

    //Destructor
    ~IdealHelmholtzLead(){};

    void to_json(rapidjson::Value &el, rapidjson::Document &doc){
        el.AddMember("type","IdealHelmholtzLead",doc.GetAllocator());
        el.AddMember("a1",a1,doc.GetAllocator());
        el.AddMember("a2",a2,doc.GetAllocator());
    };

    // Term and its derivatives
    double base(double tau, double delta){return log(delta)+a1+a2*tau;};
    double dTau(double tau, double delta){return a2;};
    double dTau2(double tau, double delta){return 0.0;};
    double dDelta(double tau, double delta){return 1.0/delta;};
    double dDelta2(double tau, double delta){return -1.0/delta/delta;};
    double dDelta2_dTau(double tau, double delta){return 0.0;};
    double dDelta_dTau(double tau, double delta){return 0.0;};
    double dDelta_dTau2(double tau, double delta){return 0.0;};
    double dTau3(double tau, double delta){return 0.0;};
    double dDelta3(double tau, double delta){return 2/delta/delta/delta;};
};

/// The term in the EOS used to shift the reference state of the fluid
/**
\f[
\alpha_0 = a_1+a_2\tau
\f]
*/
class IdealHelmholtzEnthalpyEntropyOffset : public BaseHelmholtzTerm
{
private:
    double a1,a2; // Use these variables internally
public:
    // Constructor
    IdealHelmholtzEnthalpyEntropyOffset(double a1, double a2){a1=a1; a2=a2;};

    //Destructor
    ~IdealHelmholtzEnthalpyEntropyOffset(){};

    void to_json(rapidjson::Value &el, rapidjson::Document &doc){
        el.AddMember("type","IdealHelmholtzEnthalpyEntropyOffset",doc.GetAllocator());
        el.AddMember("a1",a1,doc.GetAllocator());
        el.AddMember("a2",a2,doc.GetAllocator());
    };

    // Term and its derivatives
    double base(double tau, double delta){return a1+a2*tau;};
    double dTau(double tau, double delta){return a2;};
    double dTau2(double tau, double delta){return 0.0;};
    double dDelta(double tau, double delta){return 0.0;};
    double dDelta2(double tau, double delta){return 0.0;};
    double dDelta2_dTau(double tau, double delta){return 0.0;};
    double dDelta_dTau(double tau, double delta){return 0.0;};
    double dDelta_dTau2(double tau, double delta){return 0.0;};
    double dTau3(double tau, double delta){return 0.0;};
    double dDelta3(double tau, double delta){return 0.0;};
};


/**
\f[
\alpha_0 = a_1\ln\tau
\f]
*/
class IdealHelmholtzLogTau : public BaseHelmholtzTerm
{
private:
    double a1;
public:
    // Constructor
    IdealHelmholtzLogTau(double a1){this->a1=a1;};

    //Destructor
    ~IdealHelmholtzLogTau(){};

    void to_json(rapidjson::Value &el, rapidjson::Document &doc){
        el.AddMember("type","IdealHelmholtzLogTau",doc.GetAllocator());
        el.AddMember("a1",a1,doc.GetAllocator());
    };

    // Term and its derivatives
    double base(double tau, double delta){return a1*log(tau);};
    double dTau(double tau, double delta){return a1/tau;};
    double dTau2(double tau, double delta){return -a1/tau/tau;};
    double dTau3(double tau, double delta){return 2*a1/tau/tau/tau;};
    double dDelta(double tau, double delta){return 0.0;};
    double dDelta2(double tau, double delta){return 0.0;};
    double dDelta2_dTau(double tau, double delta){return 0.0;};
    double dDelta_dTau(double tau, double delta){return 0.0;};
    double dDelta_dTau2(double tau, double delta){return 0.0;};
    double dDelta3(double tau, double delta){return 0.0;};
};


#endif