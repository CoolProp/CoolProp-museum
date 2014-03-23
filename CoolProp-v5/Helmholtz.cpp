#include <numeric>
#include "Helmholtz.h"

namespace CoolProp{
////std::vector<double> ResidualHelmholtzTerm::dDelta(std::vector<double> tau, std::vector<double> delta)
////{
////    std::vector<double> out = tau;
////    for (int i = 0; i < (int)tau.size(); i++)
////    {
////        out[i] = dDelta(tau[i],delta[i]);
////    }
////    return out;
////}
////std::vector<double> ResidualHelmholtzTerm::dDelta2(std::vector<double> tau, std::vector<double> delta)
////{
////    std::vector<double> out = tau;
////    for (int i = 0; i < (int)tau.size(); i++)
////    {
////        out[i] = dDelta2(tau[i],delta[i]);
////    }
////    return out;
////}
////std::vector<double> ResidualHelmholtzTerm::dTau2(std::vector<double> tau, std::vector<double> delta)
////{
////    std::vector<double> out = tau;
////    for (int i = 0; i < (int)tau.size(); i++)
////    {
////        out[i] = dTau2(tau[i],delta[i]);
////    }
////    return out;
////}
////std::vector<double> ResidualHelmholtzTerm::dDelta_dTau(std::vector<double> tau, std::vector<double> delta)
////{
////    std::vector<double> out = tau;
////    for (int i = 0; i < (int)tau.size(); i++)
////    {
////        out[i] = dDelta_dTau(tau[i],delta[i]);
////    }
////    return out;
////}
//

void ResidualHelmholtzPower::to_json(rapidjson::Value &el, rapidjson::Document &doc)
{
    el.AddMember("type","ResidualHelmholtzPower",doc.GetAllocator());
    rapidjson::Value _n(rapidjson::kArrayType), _d(rapidjson::kArrayType), _t(rapidjson::kArrayType), _l(rapidjson::kArrayType);
    for (unsigned int i = 0; i <= N; ++i)
    {
        ResidualHelmholtzPowerElement &el = elements[i];
        _n.PushBack((double)el.n,doc.GetAllocator());
        _d.PushBack((double)el.d,doc.GetAllocator());
        _t.PushBack((double)el.t,doc.GetAllocator());
        _l.PushBack((double)el.l,doc.GetAllocator());
    }
    el.AddMember("n",_n,doc.GetAllocator());
    el.AddMember("d",_d,doc.GetAllocator());
    el.AddMember("t",_t,doc.GetAllocator());
    el.AddMember("l",_l,doc.GetAllocator());
}
long double ResidualHelmholtzPower::base(const long double &tau, const long double &delta)
{
    double log_tau = log(tau), log_delta = log(delta);
    for (std::size_t i = 0; i < N; ++i)
    {
        ResidualHelmholtzPowerElement &el = elements[i];
        long double ni = el.n, di = el.d, ti = el.t, lid = el.ld;
        int li = el.l;
        if (li > 0){
            s[i] = ni*exp(ti*log_tau+di*log_delta-pow(delta,li));
        }
        else{
            s[i] = ni*exp(ti*log_tau+di*log_delta);
        }
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzPower::dDelta(const long double &tau, const long double &delta)
{
    double log_tau = log(tau), log_delta = log(delta);
    for (std::size_t i = 0; i < N; ++i)
    {
        ResidualHelmholtzPowerElement &el = elements[i];
        long double ni = el.n, di = el.d, ti = el.t, lid = el.ld, pow_delta_li;
        int li = el.l;
        if (li > 0){
            pow_delta_li = pow(delta, li);
            s[i] = ni*(di-lid*pow_delta_li)*exp(ti*log_tau+(di-1)*log_delta-pow_delta_li);
        }
        else{
            s[i] = ni*di*exp(ti*log_tau+(di-1)*log_delta);
        }
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzPower::dTau(const long double &tau, const long double &delta)
{
    double log_tau = log(tau), log_delta = log(delta);
    for (std::size_t i = 0; i < N; ++i)
    {
        ResidualHelmholtzPowerElement &el = elements[i];
        long double ni = el.n, di = el.d, ti = el.t, lid = el.ld;
        int li = el.l;
        if (li>0)
        {
            s[i] = ni*ti*exp((ti-1)*log_tau+di*log_delta-pow(delta, li));
        }
        else
            s[i] = ni*ti*exp((ti-1)*log_tau+di*log_delta);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};

long double ResidualHelmholtzPower::dDelta2(const long double &tau, const long double &delta)
{
    double log_tau = log(tau), log_delta = log(delta), pow_delta_li;
    for (std::size_t i = 0; i < N; ++i)
    {
        ResidualHelmholtzPowerElement &el = elements[i];
        long double ni = el.n, di = el.d, ti = el.t, lid = el.ld;
        int li = el.l;
        if (li>0)
        {
            pow_delta_li = pow(delta, li);
            s[i] = ni*((di-li*pow_delta_li)*(di-1.0-li*pow_delta_li) - li*li*pow_delta_li)*exp(ti*log_tau+(di-2)*log_delta-pow_delta_li);
        }
        else
            s[i] = ni*di*(di-1.0)*exp(ti*log_tau+(di-2)*log_delta);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzPower::dDelta_dTau(const long double &tau, const long double &delta)
{
    double log_tau = log(tau), log_delta = log(delta), pow_delta_li;
    for (std::size_t i = 0; i < N; ++i)
    {
        ResidualHelmholtzPowerElement &el = elements[i];
        long double ni = el.n, di = el.d, ti = el.t, lid = el.ld;
        int li = el.l;
        if (li>0)
        {
            pow_delta_li = pow(delta, li);
            s[i] = ni*ti*(di-li*pow_delta_li)*exp((ti-1)*log_tau+(di-1)*log_delta-pow_delta_li);
        }
        else
            s[i] = ni*di*ti*exp((ti-1)*log_tau+(di-1)*log_delta);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzPower::dTau2(const long double &tau, const long double &delta)
{
    double log_tau = log(tau), log_delta = log(delta), pow_delta_li;
    for (std::size_t i = 0; i < N; ++i)
    {
        ResidualHelmholtzPowerElement &el = elements[i];
        long double ni = el.n, di = el.d, ti = el.t, lid = el.ld;
        int li = el.l;
        if (li>0)
        {
            pow_delta_li = pow(delta, li);
            s[i] = ni*ti*(ti-1)*exp((ti-2)*log_tau+di*log_delta-pow_delta_li);
        }
        else
            s[i] = ni*ti*(ti-1)*exp((ti-2)*log_tau+di*log_delta);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};


long double ResidualHelmholtzPower::dDelta3(const long double &tau, const long double &delta)
{
    double log_tau = log(tau), log_delta = log(delta), pow_delta_li;
    for (std::size_t i = 0; i < N; ++i)
    {
        ResidualHelmholtzPowerElement &el = elements[i];
        long double ni = el.n, di = el.d, ti = el.t, lid = el.ld;
        int li = el.l;
        if (li>0)
        {
            pow_delta_li = pow(delta, li);
            double bracket = (di*(di-1)*(di-2)+pow_delta_li*(-2*li+6*di*li-3*di*di*li-3*di*li*li+3*li*li-li*li*li)+pow_delta_li*pow_delta_li*(3*di*li*li-3*li*li+3*li*li*li)-li*li*li*pow_delta_li*pow_delta_li*pow_delta_li);
            s[i] = ni*bracket*exp(ti*log_tau+(di-3)*log_delta-pow_delta_li);
        }
        else
            s[i] = ni*di*(di-1.0)*(di-2)*exp(ti*log_tau+(di-3)*log_delta);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzPower::dDelta2_dTau(const long double &tau, const long double &delta)
{
    double log_tau = log(tau), log_delta = log(delta), pow_delta_li;
    for (std::size_t i = 0; i < N; ++i)
    {
        ResidualHelmholtzPowerElement &el = elements[i];
        long double ni = el.n, di = el.d, ti = el.t, lid = el.ld;
        int li = el.l;
        if (li>0)
        {
            pow_delta_li = pow(delta, li);
            s[i] = ni*ti*(((di-li*pow_delta_li))*(di-1-li*pow_delta_li)-li*li*pow_delta_li)*exp((ti-1)*log_tau+(di-2)*log_delta-pow_delta_li);
        }
        else
            s[i] = ni*di*ti*(di-1)*exp((ti-1)*log_tau+(di-2)*log_delta);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzPower::dDelta_dTau2(const long double &tau, const long double &delta)
{
    double log_tau = log(tau), log_delta = log(delta), pow_delta_li;
    for (std::size_t i = 0; i < N; ++i)
    {
        ResidualHelmholtzPowerElement &el = elements[i];
        long double ni = el.n, di = el.d, ti = el.t, lid = el.ld;
        int li = el.l;
        if (li>0)
        {
            pow_delta_li = pow(delta, li);
            s[i] = ni*ti*(ti-1)*(di-li*pow_delta_li)*exp((ti-2)*log_tau+(di-1)*log_delta-pow_delta_li);
        }
        else
            s[i] = ni*ti*(ti-1)*di*exp((ti-2)*log_tau+(di-1)*log_delta);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzPower::dTau3(const long double &tau, const long double &delta)
{
    double log_tau = log(tau), log_delta = log(delta);
    for (std::size_t i = 0; i < N; ++i)
    {
        ResidualHelmholtzPowerElement &el = elements[i];
        long double ni = el.n, di = el.d, ti = el.t, lid = el.ld;
        int li = el.l;
        if (li > 0)
        {
            s[i] = ni*ti*(ti-1)*(ti-2)*exp((ti-3)*log_tau+di*log_delta-pow(delta,li));
        }
        else
            s[i] = ni*ti*(ti-1)*(ti-2)*exp((ti-3)*log_tau+di*log_delta);
    }
    return std::accumulate(s.begin(), s.end(), 0.0);
};

////void ResidualHelmholtzExponential::to_json(rapidjson::Value &el, rapidjson::Document &doc)
////{
////    el.AddMember("type","alphar_exponential",doc.GetAllocator());
////    rapidjson::Value _n(rapidjson::kArrayType), _d(rapidjson::kArrayType), _t(rapidjson::kArrayType), _l(rapidjson::kArrayType), _g(rapidjson::kArrayType);
////    for (unsigned int i=0; i<=n.size(); ++i)
////    {
////        _n.PushBack(n[i],doc.GetAllocator());
////        _d.PushBack(d[i],doc.GetAllocator());
////        _t.PushBack(t[i],doc.GetAllocator());
////        _l.PushBack(l[i],doc.GetAllocator());
////        _g.PushBack(g[i],doc.GetAllocator());
////    }
////    el.AddMember("n",_n,doc.GetAllocator());
////    el.AddMember("d",_d,doc.GetAllocator());
////    el.AddMember("t",_t,doc.GetAllocator());
////    el.AddMember("l",_l,doc.GetAllocator());
////    el.AddMember("g",_g,doc.GetAllocator());
////}
////
////// Term and its derivatives
////double ResidualHelmholtzExponential::A(double log_tau, double tau, double log_delta, double delta, int i)
////{
////    return exp(t[i]*log_tau+d[i]*log_delta-g[i]*pow(delta,l[i]));
////}
////double ResidualHelmholtzExponential::dA_dTau(double log_tau, double tau, double log_delta, double delta, int i)
////{
////    return t[i]*exp((t[i]-1)*log_tau+d[i]*log_delta-g[i]*pow(delta,l[i]));
////}
////double ResidualHelmholtzExponential::d2A_dTau2(double log_tau, double tau, double log_delta, double delta, int i)
////{
////    return t[i]*(t[i]-1)*exp((t[i]-2)*log_tau+d[i]*log_delta-g[i]*pow(delta,l[i]));
////}
////double ResidualHelmholtzExponential::d3A_dTau3(double log_tau, double tau, double log_delta, double delta, int i)
////{
////    return t[i]*(t[i]-1)*(t[i]-2)*exp((t[i]-3)*log_tau+d[i]*log_delta-g[i]*pow(delta,l[i]));
////}
////double ResidualHelmholtzExponential::d3A_dDelta_dTau2(double log_tau, double tau, double log_delta, double delta, int i)
////{
////    double pow_delta_li = pow(delta,l[i]);
////    return t[i]*(t[i]-1)*(d[i]-g[i]*l[i]*pow_delta_li)*exp((t[i]-2)*log_tau+(d[i]-1)*log_delta-g[i]*pow_delta_li);
////}
////double ResidualHelmholtzExponential::dA_dDelta(double log_tau, double tau, double log_delta, double delta, int i)
////{
////    double pow_delta_li = pow(delta,l[i]);
////    return (d[i]-g[i]*l[i]*pow_delta_li)*exp(t[i]*log_tau+(d[i]-1)*log_delta-g[i]*pow_delta_li);
////    
////}
////double ResidualHelmholtzExponential::d2A_dDelta2(double log_tau, double tau, double log_delta, double delta, int i)
////{
////    double pow_delta_li = pow(delta,l[i]);
////    // Typo in Span, 2000, re-derived from Sympy
////    double bracket = d[i]*d[i] - 2*d[i]*pow(delta,l[i])*g[i]*l[i] - d[i] + pow(delta,2*l[i])*g[i]*g[i]*l[i]*l[i] - pow(delta,l[i])*g[i]*l[i]*l[i] + pow(delta,l[i])*g[i]*l[i];
////    return bracket*exp(t[i]*log_tau+(d[i]-2)*log_delta-g[i]*pow_delta_li);
////}
////double ResidualHelmholtzExponential::d3A_dDelta3(double log_tau, double tau, double log_delta, double delta, int i)
////{
////    // >> n_i, tau, t_i, d_i, delta, g_i, l_i = symbols(' n_i tau t_i d_i delta g_i l_i')
////    // >> phir = n_i*tau**t_i*delta**d_i*exp(-g_i*pow(delta,l_i))
////    // >> simplify(diff(diff(diff(phir,delta),delta),delta))
////    double pow_delta_li = pow(delta,l[i]);
////    double pow_delta_2li = pow(delta,2*l[i]);
////    double pow_delta_3li = pow(delta,3*l[i]);
////    double bracket = d[i]*d[i]*d[i] - 3*d[i]*d[i]*pow_delta_li*g[i]*l[i] - 3*d[i]*d[i] + 3*d[i]*pow_delta_2li*g[i]*g[i]*l[i]*l[i] - 3*d[i]*pow_delta_li*g[i]*l[i]*l[i] + 6*d[i]*pow_delta_li*g[i]*l[i] + 2*d[i] - pow_delta_3li*g[i]*g[i]*g[i]*l[i]*l[i]*l[i] + 3*pow_delta_2li*g[i]*g[i]*l[i]*l[i]*l[i] - 3*pow_delta_2li*g[i]*g[i]*l[i]*l[i] - pow_delta_li*g[i]*l[i]*l[i]*l[i] + 3*pow_delta_li*g[i]*l[i]*l[i] - 2*pow_delta_li*g[i]*l[i];
////    return bracket*exp(t[i]*log_tau+(d[i]-3)*log_delta-g[i]*pow_delta_li);
////}
////double ResidualHelmholtzExponential::d3A_dDelta2_dTau(double log_tau, double tau, double log_delta, double delta, int i)
////{
////    double pow_delta_li = pow(delta,l[i]);
////    // Typo in Span, 2000, re-derived from Sympy
////    double bracket = d[i]*d[i] - 2*d[i]*pow(delta,l[i])*g[i]*l[i] - d[i] + pow(delta,2*l[i])*g[i]*g[i]*l[i]*l[i] - pow(delta,l[i])*g[i]*l[i]*l[i] + pow(delta,l[i])*g[i]*l[i];
////    return t[i]*bracket*exp((t[i]-1)*log_tau+(d[i]-2)*log_delta-g[i]*pow_delta_li);
////}
////double ResidualHelmholtzExponential::d2A_dDelta_dTau(double log_tau, double tau, double log_delta, double delta, int i)
////{
////    double pow_delta_li = pow(delta,l[i]);
////    return t[i]*(d[i]-g[i]*l[i]*pow_delta_li)*exp((t[i]-1)*log_tau+(d[i]-1)*log_delta-g[i]*pow_delta_li);
////}


void ResidualHelmholtzGaussian::to_json(rapidjson::Value &el, rapidjson::Document &doc)
{
    el.AddMember("type","ResidualHelmholtzGaussian",doc.GetAllocator());
    rapidjson::Value _n(rapidjson::kArrayType), _d(rapidjson::kArrayType), _t(rapidjson::kArrayType), 
        _eta(rapidjson::kArrayType), _epsilon(rapidjson::kArrayType), _beta(rapidjson::kArrayType), _gamma(rapidjson::kArrayType);
    for (unsigned int i=0; i<=N; ++i)
    {
        ResidualHelmholtzGaussianElement &el = elements[i];
        _n.PushBack((double)el.n,doc.GetAllocator());
        _d.PushBack((double)el.d,doc.GetAllocator());
        _t.PushBack((double)el.t,doc.GetAllocator());
        _eta.PushBack((double)el.eta,doc.GetAllocator());
        _epsilon.PushBack((double)el.epsilon,doc.GetAllocator());
        _beta.PushBack((double)el.beta,doc.GetAllocator());
        _gamma.PushBack((double)el.gamma,doc.GetAllocator());
    }
    el.AddMember("n",_n,doc.GetAllocator());
    el.AddMember("d",_d,doc.GetAllocator());
    el.AddMember("t",_t,doc.GetAllocator());
    el.AddMember("eta",_eta,doc.GetAllocator());
    el.AddMember("epsilon",_epsilon,doc.GetAllocator());
    el.AddMember("beta",_beta,doc.GetAllocator());
    el.AddMember("gamma",_gamma,doc.GetAllocator());
}
long double ResidualHelmholtzGaussian::base(const long double &tau, const long double &delta)
{
	double log_tau = log(tau), log_delta = log(delta);
	for (std::size_t i=0; i<N; ++i)
	{
        ResidualHelmholtzGaussianElement &el = elements[i];
        double psi=exp(-el.eta*pow(delta-el.epsilon,2)-el.beta*pow(tau-el.gamma,2));
        s[i] = el.n*pow(delta,el.d)*pow(tau,el.t)*psi;
	}
	return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzGaussian::dDelta(const long double &tau, const long double &delta)
{
	double log_tau = log(tau), log_delta = log(delta);
	for (std::size_t i=0; i<N; ++i)
	{
        ResidualHelmholtzGaussianElement &el = elements[i];
        double psi=exp(-el.eta*pow(delta-el.epsilon,2)-el.beta*pow(tau-el.gamma,2));
        s[i] = el.n*pow(delta,el.d)*pow(tau,el.t)*psi*(el.d/delta-2.0*el.eta*(delta-el.epsilon));
	}
	return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzGaussian::dTau(const long double &tau, const long double &delta)
{
	double log_tau = log(tau), log_delta = log(delta);
	for (std::size_t i=0; i<N; ++i)
	{
        ResidualHelmholtzGaussianElement &el = elements[i];
        double psi=exp(-el.eta*pow(delta-el.epsilon,2)-el.beta*pow(tau-el.gamma,2));
        s[i] = el.n*pow(delta,el.d)*pow(tau,el.t)*psi*(el.t/tau-2.0*el.beta*(tau-el.gamma));
	}
	return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzGaussian::dDelta2(const long double &tau, const long double &delta)
{
	double log_tau = log(tau), log_delta = log(delta);
	for (std::size_t i=0; i<N; ++i)
	{
        ResidualHelmholtzGaussianElement &el = elements[i];
        double psi=exp(-el.eta*pow(delta-el.epsilon,2)-el.beta*pow(tau-el.gamma,2));
        s[i] = el.n*pow(tau,el.t)*psi*(-2.0*el.eta*pow(delta,el.d)+4.0*pow(el.eta,2)*pow(delta,el.d)*pow(delta-el.epsilon,2)-4.0*el.d*el.eta*pow(delta,el.d-1)*(delta-el.epsilon)+el.d*(el.d-1.0)*pow(delta,el.d-2));
	}
	return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzGaussian::dDelta_dTau(const long double &tau, const long double &delta)
{
	double log_tau = log(tau), log_delta = log(delta);
	for (std::size_t i=0; i<N; ++i)
	{
        ResidualHelmholtzGaussianElement &el = elements[i];
        double psi=exp(-el.eta*pow(delta-el.epsilon,2)-el.beta*pow(tau-el.gamma,2));
        s[i] = el.n*pow(delta,el.d)*pow(tau,el.t)*psi*(el.d/delta-2.0*el.eta*(delta-el.epsilon))*(el.t/tau-2.0*el.beta*(tau-el.gamma));
	}
	return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzGaussian::dTau2(const long double &tau, const long double &delta)
{
	double log_tau = log(tau), log_delta = log(delta);
	for (std::size_t i=0; i<N; ++i)
	{
        ResidualHelmholtzGaussianElement &el = elements[i];
        double psi=exp(-el.eta*pow(delta-el.epsilon,2)-el.beta*pow(tau-el.gamma,2));
        s[i] = el.n*pow(delta,el.d)*pow(tau,el.t)*psi*(pow(el.t/tau-2.0*el.beta*(tau-el.gamma),2)-el.t/pow(tau,2)-2.0*el.beta);
	}
	return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzGaussian::dDelta3(const long double &tau, const long double &delta)
{
	double log_tau = log(tau), log_delta = log(delta);
	for (std::size_t i=0; i<N; ++i)
	{
        ResidualHelmholtzGaussianElement &el = elements[i];
        double psi=exp(-el.eta*pow(delta-el.epsilon,2)-el.beta*pow(tau-el.gamma,2));
        double bracket = pow(el.t/tau-2.0*el.beta*(tau-el.gamma),2)-el.t/pow(tau,2)-2.0*el.beta;
        s[i] = el.n*pow(tau,el.t)*pow(delta,el.d-3)*psi*bracket;
	}
	return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzGaussian::dDelta2_dTau(const long double &tau, const long double &delta)
{
	double log_tau = log(tau), log_delta = log(delta);
	for (std::size_t i=0; i<N; ++i)
	{
        ResidualHelmholtzGaussianElement &el = elements[i];
        double psi=exp(-el.eta*pow(delta-el.epsilon,2)-el.beta*pow(tau-el.gamma,2));
        s[i] = el.n*pow(tau,el.t)*psi*(el.t/tau-2.0*el.beta*(tau-el.gamma))*(-2.0*el.eta*pow(delta,el.d)+4.0*pow(el.eta,2)*pow(delta,el.d)*pow(delta-el.epsilon,2)-4.0*el.d*el.eta*pow(delta,el.d-1)*(delta-el.epsilon)+el.d*(el.d-1.0)*pow(delta,el.d-2));
	}
	return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzGaussian::dDelta_dTau2(const long double &tau, const long double &delta)
{
	double log_tau = log(tau), log_delta = log(delta);
	for (std::size_t i=0; i<N; ++i)
	{
        ResidualHelmholtzGaussianElement &el = elements[i];
        double psi=exp(-el.eta*pow(delta-el.epsilon,2)-el.beta*pow(tau-el.gamma,2));
        s[i] = el.n*pow(delta,el.d)*pow(tau,el.t)*psi*(el.d/delta-2.0*el.eta*(delta-el.epsilon))*(pow(el.t-2.0*el.beta*tau*(tau-el.gamma),2)-el.t-2*el.beta*tau*tau)/tau/tau;
	}
	return std::accumulate(s.begin(), s.end(), 0.0);
};
long double ResidualHelmholtzGaussian::dTau3(const long double &tau, const long double &delta)
{
	double log_tau = log(tau), log_delta = log(delta);
	for (std::size_t i=0; i<N; ++i)
	{
        ResidualHelmholtzGaussianElement &el = elements[i];
        double psi=exp(-el.eta*pow(delta-el.epsilon,2)-el.beta*pow(tau-el.gamma,2));
        double dpsi_dTau = exp(-el.eta*pow(delta-el.epsilon,2)-el.beta*pow(tau-el.gamma,2))*(-2*el.beta*(tau-el.gamma));

        double bracket = pow(el.t/tau-2.0*el.beta*(tau-el.gamma),2)-el.t/pow(tau,2)-2.0*el.beta;
        double dbracket_dTau = 2*(el.t/tau-2.0*el.beta*(tau-el.gamma))*(-el.t/tau/tau-2*el.beta)+2*el.t/pow(tau,3);
        s[i] = el.n*pow(delta,el.d)*(el.t*pow(tau,el.t-1)*psi*bracket+pow(tau,el.t)*dpsi_dTau*bracket+pow(tau,el.t)*psi*dbracket_dTau);
	}
	return std::accumulate(s.begin(), s.end(), 0.0);
};

////void ResidualHelmholtzLemmon2005::to_json(rapidjson::Value &el, rapidjson::Document &doc)
////{
////    el.AddMember("type","alphar_Lemmon2005",doc.GetAllocator());
////    rapidjson::Value _n(rapidjson::kArrayType), _d(rapidjson::kArrayType), _t(rapidjson::kArrayType), _l(rapidjson::kArrayType), _m(rapidjson::kArrayType);
////    for (unsigned int i=0;i<=n.size();i++)
////    {
////        _n.PushBack(n[i],doc.GetAllocator());
////        _d.PushBack(d[i],doc.GetAllocator());
////        _t.PushBack(t[i],doc.GetAllocator());
////        _l.PushBack(l[i],doc.GetAllocator());
////        _m.PushBack(m[i],doc.GetAllocator());
////    }
////    el.AddMember("n",_n,doc.GetAllocator());
////    el.AddMember("d",_d,doc.GetAllocator());
////    el.AddMember("t",_t,doc.GetAllocator());
////    el.AddMember("l",_l,doc.GetAllocator());
////    el.AddMember("m",_m,doc.GetAllocator());
////}
////
////// Term and its derivatives
////double ResidualHelmholtzLemmon2005::A(double log_tau, double tau, double log_delta, double delta, int i)
////{
////    if (l[i] != 0 && m[i] != 0)
////        return exp(t[i]*log_tau+d[i]*log_delta-pow(delta,l[i])-pow(tau,m[i]));
////    else if (l[i] != 0 && m[i] == 0)
////        return exp(t[i]*log_tau+d[i]*log_delta-pow(delta,l[i]));
////    else
////        return exp(t[i]*log_tau+d[i]*log_delta);
////}
////double ResidualHelmholtzLemmon2005::dA_dTau(double log_tau, double tau, double log_delta, double delta, int i)
////{
////    if (l[i] != 0 && m[i] != 0){
////        double pow_tau_mi = pow(tau,m[i]);
////        return (t[i]-m[i]*pow_tau_mi)*exp((t[i]-1)*log_tau+d[i]*log_delta-pow(delta,l[i])-pow_tau_mi);
////    }
////    else if (l[i] != 0 && m[i] == 0)
////        return t[i]*exp((t[i]-1)*log_tau+d[i]*log_delta-pow(delta,l[i]));
////    else
////        return t[i]*exp((t[i]-1)*log_tau+d[i]*log_delta);
////}
////double ResidualHelmholtzLemmon2005::d2A_dTau2(double log_tau, double tau, double log_delta, double delta, int i)
////{
////    if (l[i] != 0 && m[i] != 0){
////        double pow_tau_mi = pow(tau,m[i]);
////        double bracket = (t[i]-m[i]*pow_tau_mi)*(t[i]-1-m[i]*pow_tau_mi)-m[i]*m[i]*pow_tau_mi;
////        return bracket*exp((t[i]-2)*log_tau+d[i]*log_delta-pow(delta,l[i])-pow_tau_mi);
////    }
////    else if (l[i] != 0 && m[i] == 0)
////        return t[i]*(t[i]-1)*exp((t[i]-2)*log_tau+d[i]*log_delta-pow(delta,l[i]));
////    else
////        return t[i]*(t[i]-1)*exp((t[i]-2)*log_tau+d[i]*log_delta);
////}
////
/////**
////
////\f[
////\frac{{{\partial ^2}{\alpha ^r}}}{{\partial {\tau ^2}}} = {N_k}{\delta ^{{d_k}}}{\tau ^{{t_k} - 2}}\exp \left( { - {\delta ^{{l_k}}}} \right)\exp \left( { - {\tau ^{{m_k}}}} \right)\left[ {\left( {{t_k} - {m_k}{\tau ^{{m_k}}}} \right)\left( {{t_k} - 1 - {m_k}{\tau ^{{m_k}}}} \right) - m_k^2{\tau ^{{m_k}}}} \right]\\
////\f]
////\f[
////\frac{{{\partial ^2}{\alpha ^r}}}{{\partial {\tau ^2}}} = {N_k}{\delta ^{{d_k}}}\exp \left( { - {\delta ^{{l_k}}}} \right){\tau ^{{t_k} - 2}}\exp \left( { - {\tau ^{{m_k}}}} \right)\left[ {\left( {{t_k} - {m_k}{\tau ^{{m_k}}}} \right)\left( {{t_k} - 1 - {m_k}{\tau ^{{m_k}}}} \right) - m_k^2{\tau ^{{m_k}}}} \right]\\
////\f]
////Group all the terms that don't depend on \$ \tau \$
////\f[
////\frac{{{\partial ^2}{\alpha ^r}}}{{\partial {\tau ^2}}} = A{\tau ^{{t_k} - 2}}\exp \left( { - {\tau ^{{m_k}}}} \right)\left[ {\left( {{t_k} - {m_k}{\tau ^{{m_k}}}} \right)\left( {{t_k} - 1 - {m_k}{\tau ^{{m_k}}}} \right) - m_k^2{\tau ^{{m_k}}}} \right]\\
////\f]
////\f[
////\frac{1}{A}\frac{{{\partial ^3}{\alpha ^r}}}{{\partial {\tau ^3}}} = {\tau ^{{t_k} - 2}}\exp \left( { - {\tau ^{{m_k}}}} \right)\frac{\partial }{{\partial \tau }}\left[ {\left( {{t_k} - {m_k}{\tau ^{{m_k}}}} \right)\left( {{t_k} - 1 - {m_k}{\tau ^{{m_k}}}} \right) - m_k^2{\tau ^{{m_k}}}} \right] + \frac{\partial }{{\partial \tau }}\left[ {{\tau ^{{t_k} - 2}}\exp \left( { - {\tau ^{{m_k}}}} \right)} \right]\left[ {\left( {{t_k} - {m_k}{\tau ^{{m_k}}}} \right)\left( {{t_k} - 1 - {m_k}{\tau ^{{m_k}}}} \right) - m_k^2{\tau ^{{m_k}}}} \right]\\
////\f]
////\f[
////\frac{\partial }{{\partial \tau }}\left[ {{\tau ^{{t_k} - 2}}\exp \left( { - {\tau ^{{m_k}}}} \right)} \right] = ({t_k} - 2){\tau ^{{t_k} - 3}}\exp \left( { - {\tau ^{{m_k}}}} \right) + {\tau ^{{t_k} - 2}}\exp \left( { - {\tau ^{{m_k}}}} \right)( - {m_k}{\tau ^{{m_k} - 1}}) = \exp \left( { - {\tau ^{{m_k}}}} \right)\left( {({t_k} - 2){\tau ^{{t_k} - 3}} - {\tau ^{{t_k} - 2}}{m_k}{\tau ^{{m_k} - 1}}} \right)\\
////\f]
////\f[
////\frac{\partial }{{\partial \tau }}\left[ {\left( {{t_k} - {m_k}{\tau ^{{m_k}}}} \right)\left( {{t_k} - 1 - {m_k}{\tau ^{{m_k}}}} \right) - m_k^2{\tau ^{{m_k}}}} \right] = \left( {{t_k} - {m_k}{\tau ^{{m_k}}}} \right)\left( { - m_k^2{\tau ^{{m_k} - 1}}} \right) + \left( { - m_k^2{\tau ^{{m_k} - 1}}} \right)\left( {{t_k} - 1 - {m_k}{\tau ^{{m_k}}}} \right) - m_k^3{\tau ^{{m_k} - 1}} =  - m_k^2{\tau ^{{m_k} - 1}}\left[ {{t_k} - {m_k}{\tau ^{{m_k}}} + {t_k} - 1 - {m_k}{\tau ^{{m_k}}} + {m_k}} \right] =  - m_k^2{\tau ^{{m_k} - 1}}\left[ {2{t_k} - 2{m_k}{\tau ^{{m_k}}} - 1 + {m_k}} \right]\\
////\f]
////\f[
////\frac{1}{A}\frac{{{\partial ^3}{\alpha ^r}}}{{\partial {\tau ^3}}} = {\tau ^{{t_k} - 2}}\exp \left( { - {\tau ^{{m_k}}}} \right)\left( { - m_k^2{\tau ^{{m_k} - 1}}\left[ {2{t_k} - 2{m_k}{\tau ^{{m_k}}} - 1 + {m_k}} \right]} \right) + \exp \left( { - {\tau ^{{m_k}}}} \right)\left( {({t_k} - 2){\tau ^{{t_k} - 3}} - {\tau ^{{t_k} - 2}}{m_k}{\tau ^{{m_k} - 1}}} \right)\left[ {\left( {{t_k} - {m_k}{\tau ^{{m_k}}}} \right)\left( {{t_k} - 1 - {m_k}{\tau ^{{m_k}}}} \right) - m_k^2{\tau ^{{m_k}}}} \right]\\
////\f]
////\f[
////\frac{1}{A}\frac{{{\partial ^3}{\alpha ^r}}}{{\partial {\tau ^3}}} = \exp \left( { - {\tau ^{{m_k}}}} \right)\left[ { - {\tau ^{{t_k} - 2}}m_k^2{\tau ^{{m_k} - 1}}\left[ {2{t_k} - 2{m_k}{\tau ^{{m_k}}} - 1 + {m_k}} \right] + \left( {({t_k} - 2){\tau ^{{t_k} - 3}} - {\tau ^{{t_k} - 2}}{m_k}{\tau ^{{m_k} - 1}}} \right)\left[ {\left( {{t_k} - {m_k}{\tau ^{{m_k}}}} \right)\left( {{t_k} - 1 - {m_k}{\tau ^{{m_k}}}} \right) - m_k^2{\tau ^{{m_k}}}} \right]} \right]
////\f]
////*/
////double ResidualHelmholtzLemmon2005::d3A_dTau3(double log_tau, double tau, double log_delta, double delta, int i)
////{
////    if (l[i] != 0 && m[i] != 0){
////        double pow_delta_li = pow(delta,l[i]);
////        double pow_tau_mi = pow(tau,m[i]);
////        double bracket = -pow(tau,t[i]+m[i]-3)*m[i]*m[i]*(2*t[i]-2*m[i]*pow_tau_mi-1-m[i])+((t[i]-2)*pow(tau,t[i]-3)-pow(tau,t[i]-2)*m[i]*pow(tau,m[i]-1))*((t[i]-m[i]*pow_tau_mi)*(t[i]-1-m[i]*pow_tau_mi)-m[i]*m[i]*pow_tau_mi);
////        return t[i]*(t[i]-1)*(t[i]-2)*exp((t[i]-3)*log_tau+d[i]*log_delta-pow_delta_li-pow_tau_mi);
////    }
////    else if (l[i] != 0 && m[i] == 0){
////        return t[i]*(t[i]-1)*(t[i]-2)*exp((t[i]-3)*log_tau+d[i]*log_delta-pow(delta,l[i]));
////    }
////    else
////        return t[i]*(t[i]-1)*(t[i]-2)*exp((t[i]-3)*log_tau+d[i]*log_delta);
////}
////double ResidualHelmholtzLemmon2005::d3A_dDelta_dTau2(double log_tau, double tau, double log_delta, double delta, int i)
////{
////    if (l[i] != 0 && m[i] != 0){
////        double pow_delta_li = pow(delta,l[i]);
////        double pow_tau_mi = pow(tau,m[i]);
////        // delta derivative of second tau derivative
////        double bracket = ((t[i]-m[i]*pow_tau_mi)*(t[i]-1-m[i]*pow_tau_mi)-m[i]*m[i]*pow_tau_mi)*(d[i]-l[i]*pow_delta_li);
////        return bracket*exp((t[i]-2)*log_tau+(d[i]-1)*log_delta-pow_delta_li-pow_tau_mi);
////    }
////    else if (l[i] != 0 && m[i] == 0){
////        double pow_delta_li = pow(delta,l[i]);
////        return t[i]*(t[i]-1)*(d[i]-l[i]*pow_delta_li)*exp((t[i]-2)*log_tau+(d[i]-1)*log_delta-pow_delta_li);
////    }
////    else
////        return t[i]*(t[i]-1)*d[i]*exp((t[i]-2)*log_tau+(d[i]-1)*log_delta);
////}
////double ResidualHelmholtzLemmon2005::dA_dDelta(double log_tau, double tau, double log_delta, double delta, int i)
////{
////    if (l[i] != 0 && m[i] != 0){
////        double pow_delta_li = pow(delta,l[i]);
////        double pow_tau_mi = pow(tau,m[i]);
////        return (d[i]-l[i]*pow_delta_li)*exp(t[i]*log_tau+(d[i]-1)*log_delta-pow_delta_li-pow_tau_mi);
////    }
////    else if (l[i]>0 && m[i] == 0){
////        double pow_delta_li = pow(delta,l[i]);
////        return (d[i]-l[i]*pow_delta_li)*exp(t[i]*log_tau+(d[i]-1)*log_delta-pow_delta_li);
////    }
////    else
////        return d[i]*exp(t[i]*log_tau+(d[i]-1)*log_delta);
////}
////double ResidualHelmholtzLemmon2005::d2A_dDelta2(double log_tau, double tau, double log_delta, double delta, int i)
////{
////    if (l[i] != 0 && m[i] != 0){	
////        double pow_delta_li = pow(delta,l[i]);
////        double pow_tau_mi = pow(tau,m[i]);
////        return ((d[i]-l[i]*pow_delta_li)*(d[i]-1.0-l[i]*pow_delta_li) - l[i]*l[i]*pow_delta_li)*exp(t[i]*log_tau+(d[i]-2)*log_delta-pow_delta_li-pow_tau_mi);
////    }
////        
////    else if (l[i] != 0 && m[i] == 0){
////        double pow_delta_li = pow(delta,l[i]);
////        return ((d[i]-l[i]*pow_delta_li)*(d[i]-1.0-l[i]*pow_delta_li) - l[i]*l[i]*pow_delta_li)*exp(t[i]*log_tau+(d[i]-2)*log_delta-pow_delta_li);
////    }
////    else
////        return d[i]*(d[i]-1.0)*exp(t[i]*log_tau+(d[i]-2)*log_delta);
////}
////double ResidualHelmholtzLemmon2005::d3A_dDelta3(double log_tau, double tau, double log_delta, double delta, int i)
////{
////    if (l[i] != 0 && m[i] != 0){
////        double pow_delta_li = pow(delta,l[i]);
////        double pow_tau_mi = pow(tau,m[i]);
////        double bracket = (d[i]*(d[i]-1)*(d[i]-2)+pow_delta_li*(-2*l[i]+6*d[i]*l[i]-3*d[i]*d[i]*l[i]-3*d[i]*l[i]*l[i]+3*l[i]*l[i]-l[i]*l[i]*l[i])+pow_delta_li*pow_delta_li*(3*d[i]*l[i]*l[i]-3*l[i]*l[i]+3*l[i]*l[i]*l[i])-l[i]*l[i]*l[i]*pow_delta_li*pow_delta_li*pow_delta_li);
////        return bracket*exp(t[i]*log_tau+(d[i]-3)*log_delta-pow_delta_li-pow_tau_mi);
////    }
////    else if (l[i] != 0 && m[i] == 0)
////    {
////        double pow_delta_li = pow(delta,l[i]);
////        double bracket = (d[i]*(d[i]-1)*(d[i]-2)+pow_delta_li*(-2*l[i]+6*d[i]*l[i]-3*d[i]*d[i]*l[i]-3*d[i]*l[i]*l[i]+3*l[i]*l[i]-l[i]*l[i]*l[i])+pow_delta_li*pow_delta_li*(3*d[i]*l[i]*l[i]-3*l[i]*l[i]+3*l[i]*l[i]*l[i])-l[i]*l[i]*l[i]*pow_delta_li*pow_delta_li*pow_delta_li);
////        return bracket*exp(t[i]*log_tau+(d[i]-3)*log_delta-pow_delta_li);
////    }
////    else
////        return d[i]*(d[i]-1.0)*(d[i]-2)*exp(t[i]*log_tau+(d[i]-3)*log_delta);
////}
////
////double ResidualHelmholtzLemmon2005::d3A_dDelta2_dTau(double log_tau, double tau, double log_delta, double delta, int i)
////{
////    if (l[i] != 0 && m[i] != 0){
////        double pow_delta_li = pow(delta,l[i]);
////        double pow_tau_mi = pow(tau,m[i]);
////        double bracket = (t[i]-m[i]*pow_tau_mi)*(((d[i]-l[i]*pow_delta_li))*(d[i]-1-l[i]*pow_delta_li)-l[i]*l[i]*pow_delta_li);
////        return bracket*exp((t[i]-1)*log_tau+(d[i]-2)*log_delta-pow_delta_li-pow_tau_mi);
////    }
////    else if (l[i] != 0 && m[i] == 0){
////        double pow_delta_li = pow(delta,l[i]);
////        double bracket = t[i]*(((d[i]-l[i]*pow_delta_li))*(d[i]-1-l[i]*pow_delta_li)-l[i]*l[i]*pow_delta_li);
////        return bracket*exp((t[i]-1)*log_tau+(d[i]-2)*log_delta-pow_delta_li);
////    }
////    else
////        return d[i]*t[i]*(d[i]-1)*exp((t[i]-1)*log_tau+(d[i]-2)*log_delta);
////}
////double ResidualHelmholtzLemmon2005::d2A_dDelta_dTau(double log_tau, double tau, double log_delta, double delta, int i)
////{
////    if (l[i] != 0 && m[i] != 0){
////        double pow_delta_li = pow(delta,l[i]);
////        double pow_tau_mi = pow(tau,m[i]);
////        return (d[i]-l[i]*pow_delta_li)*(t[i]-m[i]*pow_tau_mi)*exp((t[i]-1)*log_tau+(d[i]-1)*log_delta-pow_delta_li-pow_tau_mi);
////    }
////    else if (l[i] != 0 && m[i] == 0){
////        double pow_delta_li = pow(delta,l[i]);
////        return t[i]*(d[i]-l[i]*pow_delta_li)*exp((t[i]-1)*log_tau+(d[i]-1)*log_delta-pow_delta_li);
////    }
////    else
////        return d[i]*t[i]*exp((t[i]-1)*log_tau+(d[i]-1)*log_delta);
////}
////
////
////
////
////void ResidualHelmholtzGERG2008Gaussian::to_json(rapidjson::Value &el, rapidjson::Document &doc)
////{
////    el.AddMember("type","alphar_GERG2008_gaussian",doc.GetAllocator());
////    rapidjson::Value _n(rapidjson::kArrayType), _d(rapidjson::kArrayType), _t(rapidjson::kArrayType), 
////        _eta(rapidjson::kArrayType), _epsilon(rapidjson::kArrayType), _beta(rapidjson::kArrayType), _gamma(rapidjson::kArrayType);
////    for (unsigned int i=0; i<=n.size(); ++i)
////    {
////        _n.PushBack(n[i],doc.GetAllocator());
////        _d.PushBack(d[i],doc.GetAllocator());
////        _t.PushBack(t[i],doc.GetAllocator());
////        _eta.PushBack(eta[i],doc.GetAllocator());
////        _epsilon.PushBack(epsilon[i],doc.GetAllocator());
////        _beta.PushBack(beta[i],doc.GetAllocator());
////        _gamma.PushBack(gamma[i],doc.GetAllocator());
////    }
////    el.AddMember("n",_n,doc.GetAllocator());
////    el.AddMember("d",_d,doc.GetAllocator());
////    el.AddMember("t",_t,doc.GetAllocator());
////    el.AddMember("eta",_eta,doc.GetAllocator());
////    el.AddMember("epsilon",_epsilon,doc.GetAllocator());
////    el.AddMember("beta",_beta,doc.GetAllocator());
////    el.AddMember("gamma",_gamma,doc.GetAllocator());
////}
////
////double ResidualHelmholtzGERG2008Gaussian::A(double log_tau, double tau, double log_delta, double delta, int i)
////{
////    return pow(tau,t[i])*pow(delta,d[i])*psi(tau,delta,i);
////}
////double ResidualHelmholtzGERG2008Gaussian::dA_dDelta(double log_tau, double tau, double log_delta, double delta, int i)
////{
////    return pow(tau,t[i])*pow(delta,d[i])*psi(tau,delta,i)*(d[i]/delta-2*eta[i]*(delta-epsilon[i])-beta[i]);
////}
////double ResidualHelmholtzGERG2008Gaussian::d2A_dDelta2(double log_tau, double tau, double log_delta, double delta, int i)
////{
////    return pow(tau,t[i])*pow(delta,d[i])*psi(tau,delta,i)*(pow(d[i]/delta-2*eta[i]*(delta-epsilon[i])-beta[i],2)-d[i]/delta/delta-2*eta[i]);
////}
////double ResidualHelmholtzGERG2008Gaussian::dA_dTau(double log_tau, double tau, double log_delta, double delta, int i)
////{
////    return t[i]*pow(tau,t[i]-1)*pow(delta,d[i])*psi(tau,delta,i);
////}
////double ResidualHelmholtzGERG2008Gaussian::d2A_dTau2(double log_tau, double tau, double log_delta, double delta, int i)
////{
////    return t[i]*(t[i]-1)*pow(tau,t[i]-2)*pow(delta,d[i])*psi(tau,delta,i);
////}
////double ResidualHelmholtzGERG2008Gaussian::d2A_dDelta_dTau(double log_tau, double tau, double log_delta, double delta, int i)
////{
////    return t[i]*pow(tau,t[i]-1)*pow(delta,d[i])*psi(tau,delta,i)*(d[i]/delta-2*eta[i]*(delta-epsilon[i])-beta[i]);
////}
//
//void ResidualHelmholtzNonAnalytic::to_json(rapidjson::Value &el, rapidjson::Document &doc)
//{
//    el.AddMember("type","alphar_critical",doc.GetAllocator());
//
//    rapidjson::Value _n(rapidjson::kArrayType), _a(rapidjson::kArrayType), _b(rapidjson::kArrayType), 
//        _beta(rapidjson::kArrayType), __A(rapidjson::kArrayType), _B(rapidjson::kArrayType), _C(rapidjson::kArrayType), _D(rapidjson::kArrayType);
//    for (unsigned int i=0; i<=n.size(); ++i)
//    {
//        _n.PushBack(n[i],doc.GetAllocator());
//        _a.PushBack(a[i],doc.GetAllocator());
//        _b.PushBack(b[i],doc.GetAllocator());
//        _beta.PushBack(beta[i],doc.GetAllocator());
//        __A.PushBack(_A[i],doc.GetAllocator());
//        _B.PushBack(B[i],doc.GetAllocator());
//        _C.PushBack(C[i],doc.GetAllocator());
//        _D.PushBack(D[i],doc.GetAllocator());
//    }
//    el.AddMember("n",_n,doc.GetAllocator());
//    el.AddMember("a",_a,doc.GetAllocator());
//    el.AddMember("b",_b,doc.GetAllocator());
//    el.AddMember("beta",_beta,doc.GetAllocator());
//    el.AddMember("A",__A,doc.GetAllocator());
//    el.AddMember("B",_B,doc.GetAllocator());
//    el.AddMember("C",_C,doc.GetAllocator());
//    el.AddMember("D",_D,doc.GetAllocator());
//}
//

long double ResidualHelmholtzNonAnalytic::base(const long double &tau, const long double &delta)
{
	long double log_tau = log(tau), log_delta = log(delta);
	for (unsigned int i=0; i<N; ++i)
	{
        ResidualHelmholtzNonAnalyticElement &el = elements[i];
        long double ni = el.n, ai = el.a, bi = el.b, betai = el.beta;
        long double Ai = el.A, Bi = el.B, Ci = el.C, Di = el.D;
        long double theta=(1.0-tau)+Ai*pow(pow(delta-1.0,2),1.0/(2.0*betai));
        long double DELTA=pow(theta,2)+Bi*pow(pow(delta-1.0,2),ai);
        long double PSI=exp(-Ci*pow(delta-1.0,2)-Di*pow(tau-1.0,2));
        long double dPSI_dDelta=-2.0*Ci*(delta-1.0)*PSI;
        long double dDELTA_dDelta=(delta-1.0)*(Ai*theta*2.0/betai*pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0)+2.0*Bi*ai*pow(pow(delta-1.0,2),ai-1.0));
        
        s[i] = ni*pow(DELTA, bi)*delta*PSI;;
	}
	return std::accumulate(s.begin(), s.end(), 0.0);
}
long double ResidualHelmholtzNonAnalytic::dDelta(const long double &tau, const long double &delta)
{
	long double log_tau = log(tau), log_delta = log(delta);
	for (unsigned int i=0; i<N; ++i)
	{
        ResidualHelmholtzNonAnalyticElement &el = elements[i];
        long double ni = el.n, ai = el.a, bi = el.b, betai = el.beta;
        long double Ai = el.A, Bi = el.B, Ci = el.C, Di = el.D;
        long double dDELTAbi_dDelta;
        long double theta=(1.0-tau)+Ai*pow(pow(delta-1.0,2),1.0/(2.0*betai));
        long double DELTA=pow(theta,2)+Bi*pow(pow(delta-1.0,2),ai);
        long double PSI=exp(-Ci*pow(delta-1.0,2)-Di*pow(tau-1.0,2));
        long double dPSI_dDelta=-2.0*Ci*(delta-1.0)*PSI;
        long double dDELTA_dDelta=(delta-1.0)*(Ai*theta*2.0/betai*pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0)+2.0*Bi*ai*pow(pow(delta-1.0,2),ai-1.0));
        
        // At critical point, DELTA is 0, and 1/0^n is undefined
        if (fabs(DELTA) < 10*DBL_EPSILON)
        {
            dDELTAbi_dDelta = 0;
        }
        else{
            dDELTAbi_dDelta=bi*pow(DELTA,bi-1.0)*dDELTA_dDelta;
        }
        s[i] = ni*(pow(DELTA,bi)*(PSI+delta*dPSI_dDelta)+dDELTAbi_dDelta*delta*PSI);
	}
	return std::accumulate(s.begin(), s.end(), 0.0);
}
long double ResidualHelmholtzNonAnalytic::dTau(const long double &tau, const long double &delta)
{
	long double log_tau = log(tau), log_delta = log(delta);
	for (unsigned int i=0; i<N; ++i)
	{
        ResidualHelmholtzNonAnalyticElement &el = elements[i];
        long double ni = el.n, ai = el.a, bi = el.b, betai = el.beta;
        long double Ai = el.A, Bi = el.B, Ci = el.C, Di = el.D;
        long double dDELTAbi_dDelta;
        long double theta=(1.0-tau)+Ai*pow(pow(delta-1.0,2),1.0/(2.0*betai));
        long double DELTA=pow(theta,2)+Bi*pow(pow(delta-1.0,2),ai);
        long double PSI=exp(-Ci*pow(delta-1.0,2)-Di*pow(tau-1.0,2));
        long double dPSI_dTau=-2.0*Di*(tau-1.0)*PSI;
        long double dDELTA_dDelta=(delta-1.0)*(Ai*theta*2.0/betai*pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0)+2.0*Bi*ai*pow(pow(delta-1.0,2),ai-1.0));
        long double dDELTAbi_dTau=-2.0*theta*bi*pow(DELTA,bi-1.0);
        
        // At critical point, DELTA is 0, and 1/0^n is undefined
        if (fabs(DELTA) < 10*DBL_EPSILON)
        {
            dDELTAbi_dDelta = 0;
        }
        else{
            dDELTAbi_dDelta=bi*pow(DELTA,bi-1.0)*dDELTA_dDelta;
        }
        s[i] = ni*delta*(dDELTAbi_dTau*PSI+pow(DELTA,bi)*dPSI_dTau);
	}
	return std::accumulate(s.begin(), s.end(), 0.0);
}

long double ResidualHelmholtzNonAnalytic::dDelta2(const long double &tau, const long double &delta)
{
	long double log_tau = log(tau), log_delta = log(delta);
	for (unsigned int i=0; i<N; ++i)
	{
        ResidualHelmholtzNonAnalyticElement &el = elements[i];
        long double ni = el.n, ai = el.a, bi = el.b, betai = el.beta;
        long double Ai = el.A, Bi = el.B, Ci = el.C, Di = el.D;
        long double dDELTAbi_dDelta;
        long double theta=(1.0-tau)+Ai*pow(pow(delta-1.0,2),1.0/(2.0*betai));
        long double DELTA=pow(theta,2)+Bi*pow(pow(delta-1.0,2),ai);
        long double PSI=exp(-Ci*pow(delta-1.0,2)-Di*pow(tau-1.0,2));
        long double dPSI_dDelta=-2.0*Ci*(delta-1.0)*PSI;
        long double dDELTA_dDelta=(delta-1.0)*(Ai*theta*2.0/betai*pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0)+2.0*Bi*ai*pow(pow(delta-1.0,2),ai-1.0));

        long double dPSI2_dDelta2=(2.0*Ci*pow(delta-1.0,2)-1.0)*2.0*Ci*PSI;
        long double dDELTA2_dDelta2=1.0/(delta-1.0)*dDELTA_dDelta+pow(delta-1.0,2)*(4.0*Bi*ai*(ai-1.0)*pow(pow(delta-1.0,2),ai-2.0)+2.0*pow(Ai/betai,2)*pow(pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0),2)+Ai*theta*4.0/betai*(1.0/(2.0*betai)-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*betai)-2.0));
        long double dDELTAbi2_dDelta2=bi*(pow(DELTA,bi-1.0)*dDELTA2_dDelta2+(bi-1.0)*pow(DELTA,bi-2.0)*pow(dDELTA_dDelta,2));
        
        // At critical point, DELTA is 0, and 1/0^n is undefined
        if (fabs(DELTA) < 10*DBL_EPSILON)
        {
            dDELTAbi_dDelta = 0;
        }
        else{
            dDELTAbi_dDelta=bi*pow(DELTA,bi-1.0)*dDELTA_dDelta;
        }

        s[i] = ni*(pow(DELTA,bi)*(2.0*dPSI_dDelta+delta*dPSI2_dDelta2)+2.0*dDELTAbi_dDelta*(PSI+delta*dPSI_dDelta)+dDELTAbi2_dDelta2*delta*PSI);
	}
	return std::accumulate(s.begin(), s.end(), 0.0);
}
long double ResidualHelmholtzNonAnalytic::dDelta_dTau(const long double &tau, const long double &delta)
{
	long double log_tau = log(tau), log_delta = log(delta);
	for (unsigned int i=0; i<N; ++i)
	{
        ResidualHelmholtzNonAnalyticElement &el = elements[i];
        long double ni = el.n, ai = el.a, bi = el.b, betai = el.beta;
        long double Ai = el.A, Bi = el.B, Ci = el.C, Di = el.D;
        long double dDELTAbi_dDelta;
        long double theta=(1.0-tau)+Ai*pow(pow(delta-1.0,2),1.0/(2.0*betai));
        long double DELTA=pow(theta,2)+Bi*pow(pow(delta-1.0,2),ai);
        long double PSI=exp(-Ci*pow(delta-1.0,2)-Di*pow(tau-1.0,2));
        long double dPSI_dDelta=-2.0*Ci*(delta-1.0)*PSI;
        long double dDELTA_dDelta=(delta-1.0)*(Ai*theta*2.0/betai*pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0)+2.0*Bi*ai*pow(pow(delta-1.0,2),ai-1.0));

        long double dPSI2_dDelta2=(2.0*Ci*pow(delta-1.0,2)-1.0)*2.0*Ci*PSI;
        long double dDELTA2_dDelta2=1.0/(delta-1.0)*dDELTA_dDelta+pow(delta-1.0,2)*(4.0*Bi*ai*(ai-1.0)*pow(pow(delta-1.0,2),ai-2.0)+2.0*pow(Ai/betai,2)*pow(pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0),2)+Ai*theta*4.0/betai*(1.0/(2.0*betai)-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*betai)-2.0));
        long double dDELTAbi2_dDelta2=bi*(pow(DELTA,bi-1.0)*dDELTA2_dDelta2+(bi-1.0)*pow(DELTA,bi-2.0)*pow(dDELTA_dDelta,2));

        long double dPSI2_dDelta_dTau=4.0*Ci*Di*(delta-1.0)*(tau-1.0)*PSI;
        long double dDELTAbi2_dDelta_dTau=-Ai*bi*2.0/betai*pow(DELTA,bi-1.0)*(delta-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0)-2.0*theta*bi*(bi-1.0)*pow(DELTA,bi-2.0)*dDELTA_dDelta;
        
        long double dPSI_dTau=-2.0*Di*(tau-1.0)*PSI;
        long double dDELTAbi_dTau=-2.0*theta*bi*pow(DELTA,bi-1.0);

        // At critical point, DELTA is 0, and 1/0^n is undefined
        if (fabs(DELTA) < 10*DBL_EPSILON)
        {
            dDELTAbi_dDelta = 0;
        }
        else{
            dDELTAbi_dDelta=bi*pow(DELTA,bi-1.0)*dDELTA_dDelta;
        }

        s[i] = ni*(pow(DELTA,bi)*(dPSI_dTau+delta*dPSI2_dDelta_dTau)+delta*dDELTAbi_dDelta*dPSI_dTau+ dDELTAbi_dTau*(PSI+delta*dPSI_dDelta)+dDELTAbi2_dDelta_dTau*delta*PSI);
	}
	return std::accumulate(s.begin(), s.end(), 0.0);
}
long double ResidualHelmholtzNonAnalytic::dTau2(const long double &tau, const long double &delta)
{
	long double log_tau = log(tau), log_delta = log(delta);
	for (unsigned int i=0; i<N; ++i)
	{
        ResidualHelmholtzNonAnalyticElement &el = elements[i];
        long double ni = el.n, ai = el.a, bi = el.b, betai = el.beta;
        long double Ai = el.A, Bi = el.B, Ci = el.C, Di = el.D;
        long double theta=(1.0-tau)+Ai*pow(pow(delta-1.0,2),1.0/(2.0*betai));
        long double DELTA=pow(theta,2)+Bi*pow(pow(delta-1.0,2),ai);
        long double PSI=exp(-Ci*pow(delta-1.0,2)-Di*pow(tau-1.0,2));
        long double dPSI_dTau=-2.0*Di*(tau-1.0)*PSI;
        long double dDELTAbi_dTau=-2.0*theta*bi*pow(DELTA,bi-1.0);
        long double dPSI2_dTau2=(2.0*Di*pow(tau-1.0,2)-1.0)*2.0*Di*PSI;
        long double dDELTAbi2_dTau2=2.0*bi*pow(DELTA,bi-1.0)+4.0*pow(theta,2)*bi*(bi-1.0)*pow(DELTA,bi-2.0);

        s[i] = ni*(dDELTAbi2_dTau2*PSI+2.0*dDELTAbi_dTau*dPSI_dTau+pow(DELTA,bi)*dPSI2_dTau2);
	}
	return std::accumulate(s.begin(), s.end(), 0.0);
}

long double ResidualHelmholtzNonAnalytic::dDelta3(const long double &tau, const long double &delta)
{
	long double log_tau = log(tau), log_delta = log(delta);
	for (unsigned int i=0; i<N; ++i)
	{
        ResidualHelmholtzNonAnalyticElement &el = elements[i];
        long double ni = el.n, ai = el.a, bi = el.b, betai = el.beta;
        long double Ai = el.A, Bi = el.B, Ci = el.C, Di = el.D;
        long double dDELTAbi_dDelta;
        long double theta=(1.0-tau)+Ai*pow(pow(delta-1.0,2),1.0/(2.0*betai));
        long double DELTA=pow(theta,2)+Bi*pow(pow(delta-1.0,2),ai);
        long double PSI=exp(-Ci*pow(delta-1.0,2)-Di*pow(tau-1.0,2));
        long double dPSI_dDelta=-2.0*Ci*(delta-1.0)*PSI;
        long double dDELTA_dDelta=(delta-1.0)*(Ai*theta*2.0/betai*pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0)+2.0*Bi*ai*pow(pow(delta-1.0,2),ai-1.0));

        long double dPSI2_dDelta2=(2.0*Ci*pow(delta-1.0,2)-1.0)*2.0*Ci*PSI;
        long double dDELTA2_dDelta2=1.0/(delta-1.0)*dDELTA_dDelta+pow(delta-1.0,2)*(4.0*Bi*ai*(ai-1.0)*pow(pow(delta-1.0,2),ai-2.0)+2.0*pow(Ai/betai,2)*pow(pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0),2)+Ai*theta*4.0/betai*(1.0/(2.0*betai)-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*betai)-2.0));
        long double dDELTAbi2_dDelta2=bi*(pow(DELTA,bi-1.0)*dDELTA2_dDelta2+(bi-1.0)*pow(DELTA,bi-2.0)*pow(dDELTA_dDelta,2));
        
        long double dtheta_dDelta = Ai/(2*betai)*pow(pow(delta-1,2),1/(2*betai)-1)*2*(delta-1);
        long double dPSI3_dDelta3=2.0*Ci*PSI*(-4*Ci*Ci*pow(delta-1.0,3)+6*Ci*(delta-1));
        long double PI = 4*Bi*ai*(ai-1)*pow(pow(delta-1,2),ai-2)+2*pow(Ai/betai,2)*pow(pow(delta-1,2),1/betai-2)+4*Ai*theta/betai*(1/(2*betai)-1)*pow(pow(delta-1,2),1/(2*betai)-2);
        long double dPI_dDelta = -8*Bi*ai*(ai-1)*(ai-2)*pow(pow(delta-1,2),ai-5.0/2.0)-8*pow(Ai/betai,2)*(1/(2*betai)-1)*pow(pow(delta-1,2),1/betai-5.0/2.0)-(8*Ai*theta)/betai*(1/(2*betai)-1)*(1/(2*betai)-2)*pow(pow(delta-1,2),1/(2*betai)-5.0/2.0)+4*Ai/betai*(1/(2*betai)-1)*pow(pow(delta-1,2),1/(2*betai)-2)*dtheta_dDelta;
        long double dDELTA3_dDelta3 = 1/(delta-1)*dDELTA2_dDelta2-1/pow(delta-1,2)*dDELTA_dDelta+pow(delta-1,2)*dPI_dDelta+2*(delta-1)*PI;        
        long double dDELTAbi3_dDelta3 = bi*(pow(DELTA,bi-1)*dDELTA3_dDelta3+dDELTA2_dDelta2*(bi-1)*pow(DELTA,bi-2)*dDELTA_dDelta+(bi-1)*(pow(DELTA,bi-2)*2*dDELTA_dDelta*dDELTA2_dDelta2+pow(dDELTA_dDelta,2)*(bi-2)*pow(DELTA,bi-3)*dDELTA_dDelta));

        // At critical point, DELTA is 0, and 1/0^n is undefined
        if (fabs(DELTA) < 10*DBL_EPSILON)
        {
            dDELTAbi_dDelta = 0;
        }
        else{
            dDELTAbi_dDelta=bi*pow(DELTA,bi-1.0)*dDELTA_dDelta;
        }

        s[i] = ni*(pow(DELTA,bi)*(3.0*dPSI2_dDelta2+delta*dPSI3_dDelta3)+3.0*dDELTAbi_dDelta*(2*dPSI_dDelta+delta*dPSI2_dDelta2)+3*dDELTAbi2_dDelta2*(PSI+delta*dPSI_dDelta)+dDELTAbi3_dDelta3*PSI*delta);
	}
	return std::accumulate(s.begin(), s.end(), 0.0);
}
long double ResidualHelmholtzNonAnalytic::dDelta_dTau2(const long double &tau, const long double &delta)
{
	long double log_tau = log(tau), log_delta = log(delta);
	for (unsigned int i=0; i<N; ++i)
	{
        ResidualHelmholtzNonAnalyticElement &el = elements[i];
        long double ni = el.n, ai = el.a, bi = el.b, betai = el.beta;
        long double Ai = el.A, Bi = el.B, Ci = el.C, Di = el.D;

        long double dDELTAbi_dDelta;
        long double theta=(1.0-tau)+Ai*pow(pow(delta-1.0,2),1.0/(2.0*betai));
        long double DELTA=pow(theta,2)+Bi*pow(pow(delta-1.0,2),ai);
        long double PSI=exp(-Ci*pow(delta-1.0,2)-Di*pow(tau-1.0,2));
        long double dPSI_dDelta=-2.0*Ci*(delta-1.0)*PSI;
        long double dDELTA_dDelta=(delta-1.0)*(Ai*theta*2.0/betai*pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0)+2.0*Bi*ai*pow(pow(delta-1.0,2),ai-1.0));

        long double dDELTAbi_dTau=-2.0*theta*bi*pow(DELTA,bi-1.0);
        long double dPSI_dTau=-2.0*Di*(tau-1.0)*PSI;

        long double dPSI2_dDelta2=(2.0*Ci*pow(delta-1.0,2)-1.0)*2.0*Ci*PSI;
        long double dDELTA2_dDelta2=1.0/(delta-1.0)*dDELTA_dDelta+pow(delta-1.0,2)*(4.0*Bi*ai*(ai-1.0)*pow(pow(delta-1.0,2),ai-2.0)+2.0*pow(Ai/betai,2)*pow(pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0),2)+Ai*theta*4.0/betai*(1.0/(2.0*betai)-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*betai)-2.0));
        long double dDELTAbi2_dDelta2=bi*(pow(DELTA,bi-1.0)*dDELTA2_dDelta2+(bi-1.0)*pow(DELTA,bi-2.0)*pow(dDELTA_dDelta,2));
        
        long double dtheta_dDelta = Ai/(2*betai)*pow(pow(delta-1,2),1/(2*betai)-1)*2*(delta-1);
        long double dPSI3_dDelta3=2.0*Ci*PSI*(-4*Ci*Ci*pow(delta-1.0,3)+6*Ci*(delta-1));
        long double PI = 4*Bi*ai*(ai-1)*pow(pow(delta-1,2),ai-2)+2*pow(Ai/betai,2)*pow(pow(delta-1,2),1/betai-2)+4*Ai*theta/betai*(1/(2*betai)-1)*pow(pow(delta-1,2),1/(2*betai)-2);
        long double dPI_dDelta = -8*Bi*ai*(ai-1)*(ai-2)*pow(pow(delta-1,2),ai-5.0/2.0)-8*pow(Ai/betai,2)*(1/(2*betai)-1)*pow(pow(delta-1,2),1/betai-5.0/2.0)-(8*Ai*theta)/betai*(1/(2*betai)-1)*(1/(2*betai)-2)*pow(pow(delta-1,2),1/(2*betai)-5.0/2.0)+4*Ai/betai*(1/(2*betai)-1)*pow(pow(delta-1,2),1/(2*betai)-2)*dtheta_dDelta;
        long double dDELTA3_dDelta3 = 1/(delta-1)*dDELTA2_dDelta2-1/pow(delta-1,2)*dDELTA_dDelta+pow(delta-1,2)*dPI_dDelta+2*(delta-1)*PI;        
        long double dDELTAbi3_dDelta3 = bi*(pow(DELTA,bi-1)*dDELTA3_dDelta3+dDELTA2_dDelta2*(bi-1)*pow(DELTA,bi-2)*dDELTA_dDelta+(bi-1)*(pow(DELTA,bi-2)*2*dDELTA_dDelta*dDELTA2_dDelta2+pow(dDELTA_dDelta,2)*(bi-2)*pow(DELTA,bi-3)*dDELTA_dDelta));

        long double dPSI2_dDelta_dTau=4.0*Ci*Di*(delta-1.0)*(tau-1.0)*PSI;
        long double dDELTAbi2_dDelta_dTau=-Ai*bi*2.0/betai*pow(DELTA,bi-1.0)*(delta-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0)-2.0*theta*bi*(bi-1.0)*pow(DELTA,bi-2.0)*dDELTA_dDelta;

        long double dPSI2_dTau2=(2.0*Di*pow(tau-1.0,2)-1.0)*2.0*Di*PSI;
        long double dDELTAbi2_dTau2=2.0*bi*pow(DELTA,bi-1.0)+4.0*pow(theta,2)*bi*(bi-1.0)*pow(DELTA,bi-2.0);

        long double dPSI3_dDelta_dTau2 = 2*Di*(2*Di*pow(tau-1,2)-1)*dPSI_dDelta;
        long double dDELTAbi3_dDelta_dTau2 = 2*bi*(bi-1)*pow(DELTA,bi-2)*dDELTA_dDelta+4*pow(theta,2)*bi*(bi-1)*(bi-2)*pow(DELTA,bi-3)*dDELTA_dDelta+8*theta*bi*(bi-1)*pow(DELTA,bi-2)*dtheta_dDelta;

        // At critical point, DELTA is 0, and 1/0^n is undefined
        if (fabs(DELTA) < 10*DBL_EPSILON)
        {
            dDELTAbi_dDelta = 0;
        }
        else{
            dDELTAbi_dDelta=bi*pow(DELTA,bi-1.0)*dDELTA_dDelta;
        }

        s[i] = ni*delta*(dDELTAbi2_dTau2*dPSI_dDelta+dDELTAbi3_dDelta_dTau2*PSI+2*dDELTAbi_dTau*dPSI2_dDelta_dTau+2.0*dDELTAbi2_dDelta_dTau*dPSI_dTau+pow(DELTA,bi)*dPSI3_dDelta_dTau2+dDELTAbi_dDelta*dPSI2_dTau2)+ni*(dDELTAbi2_dTau2*PSI+2.0*dDELTAbi_dTau*dPSI_dTau+pow(DELTA,bi)*dPSI2_dTau2);
	}
	return std::accumulate(s.begin(), s.end(), 0.0);
}

long double ResidualHelmholtzNonAnalytic::dDelta2_dTau(const long double &tau, const long double &delta)
{
	long double log_tau = log(tau), log_delta = log(delta);
	for (unsigned int i=0; i<N; ++i)
	{
        ResidualHelmholtzNonAnalyticElement &el = elements[i];
        long double ni = el.n, ai = el.a, bi = el.b, betai = el.beta;
        long double Ai = el.A, Bi = el.B, Ci = el.C, Di = el.D;

        long double dDELTAbi_dDelta;
        long double theta=(1.0-tau)+Ai*pow(pow(delta-1.0,2),1.0/(2.0*betai));
        long double DELTA=pow(theta,2)+Bi*pow(pow(delta-1.0,2),ai);
        long double PSI=exp(-Ci*pow(delta-1.0,2)-Di*pow(tau-1.0,2));
        long double dPSI_dDelta=-2.0*Ci*(delta-1.0)*PSI;
        long double dDELTA_dDelta=(delta-1.0)*(Ai*theta*2.0/betai*pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0)+2.0*Bi*ai*pow(pow(delta-1.0,2),ai-1.0));

        long double dDELTAbi_dTau=-2.0*theta*bi*pow(DELTA,bi-1.0);
        long double dPSI_dTau=-2.0*Di*(tau-1.0)*PSI;

        long double dPSI2_dDelta2=(2.0*Ci*pow(delta-1.0,2)-1.0)*2.0*Ci*PSI;
        long double dDELTA2_dDelta2=1.0/(delta-1.0)*dDELTA_dDelta+pow(delta-1.0,2)*(4.0*Bi*ai*(ai-1.0)*pow(pow(delta-1.0,2),ai-2.0)+2.0*pow(Ai/betai,2)*pow(pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0),2)+Ai*theta*4.0/betai*(1.0/(2.0*betai)-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*betai)-2.0));
        long double dDELTAbi2_dDelta2=bi*(pow(DELTA,bi-1.0)*dDELTA2_dDelta2+(bi-1.0)*pow(DELTA,bi-2.0)*pow(dDELTA_dDelta,2));
        
        long double dtheta_dDelta = Ai/(2*betai)*pow(pow(delta-1,2),1/(2*betai)-1)*2*(delta-1);
        long double dPSI3_dDelta3=2.0*Ci*PSI*(-4*Ci*Ci*pow(delta-1.0,3)+6*Ci*(delta-1));
        long double PI = 4*Bi*ai*(ai-1)*pow(pow(delta-1,2),ai-2)+2*pow(Ai/betai,2)*pow(pow(delta-1,2),1/betai-2)+4*Ai*theta/betai*(1/(2*betai)-1)*pow(pow(delta-1,2),1/(2*betai)-2);
        long double dPI_dDelta = -8*Bi*ai*(ai-1)*(ai-2)*pow(pow(delta-1,2),ai-5.0/2.0)-8*pow(Ai/betai,2)*(1/(2*betai)-1)*pow(pow(delta-1,2),1/betai-5.0/2.0)-(8*Ai*theta)/betai*(1/(2*betai)-1)*(1/(2*betai)-2)*pow(pow(delta-1,2),1/(2*betai)-5.0/2.0)+4*Ai/betai*(1/(2*betai)-1)*pow(pow(delta-1,2),1/(2*betai)-2)*dtheta_dDelta;
        long double dDELTA3_dDelta3 = 1/(delta-1)*dDELTA2_dDelta2-1/pow(delta-1,2)*dDELTA_dDelta+pow(delta-1,2)*dPI_dDelta+2*(delta-1)*PI;        
        long double dDELTAbi3_dDelta3 = bi*(pow(DELTA,bi-1)*dDELTA3_dDelta3+dDELTA2_dDelta2*(bi-1)*pow(DELTA,bi-2)*dDELTA_dDelta+(bi-1)*(pow(DELTA,bi-2)*2*dDELTA_dDelta*dDELTA2_dDelta2+pow(dDELTA_dDelta,2)*(bi-2)*pow(DELTA,bi-3)*dDELTA_dDelta));

        long double dPSI2_dDelta_dTau=4.0*Ci*Di*(delta-1.0)*(tau-1.0)*PSI;
        long double dDELTAbi2_dDelta_dTau=-Ai*bi*2.0/betai*pow(DELTA,bi-1.0)*(delta-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*betai)-1.0)-2.0*theta*bi*(bi-1.0)*pow(DELTA,bi-2.0)*dDELTA_dDelta;

        long double dPSI2_dTau2=(2.0*Di*pow(tau-1.0,2)-1.0)*2.0*Di*PSI;
        long double dDELTAbi2_dTau2=2.0*bi*pow(DELTA,bi-1.0)+4.0*pow(theta,2)*bi*(bi-1.0)*pow(DELTA,bi-2.0);

        // At critical point, DELTA is 0, and 1/0^n is undefined
        if (fabs(DELTA) < 10*DBL_EPSILON)
        {
            dDELTAbi_dDelta = 0;
        }
        else{
            dDELTAbi_dDelta=bi*pow(DELTA,bi-1.0)*dDELTA_dDelta;
        }
        //Following Terms added for this derivative
        long double dPSI3_dDelta2_dTau = (2.0*Ci*pow(delta-1.0,2)-1.0)*2.0*Ci*dPSI_dTau;
        long double dDELTA_dTau = -2*theta;
        long double dDELTA2_dDelta_dTau = 2.0*Ai/(betai)*pow(pow(delta-1,2),1.0/(2.0*betai)-0.5);
        long double dDELTA3_dDelta2_dTau = 2.0*Ai*(betai-1)/(betai*betai)*pow(pow(delta-1,2),1/(2*betai)-1.0);
        
        long double dDELTAbim1_dTau = (bi-1)*pow(DELTA,bi-2)*dDELTA_dTau;
        long double dDELTAbim2_dTau = (bi-2)*pow(DELTA,bi-3)*dDELTA_dTau;
        long double Line11 = dDELTAbim1_dTau*dDELTA2_dDelta2 + pow(DELTA,bi-1)*dDELTA3_dDelta2_dTau;
        long double Line21 = (bi-1)*(dDELTAbim2_dTau*pow(dDELTA_dDelta,2)+pow(DELTA,bi-2)*2*dDELTA_dDelta*dDELTA2_dDelta_dTau);
        long double dDELTAbi3_dDelta2_dTau = bi*(Line11+Line21);

        long double Line1 = pow(DELTA,bi)*(2*dPSI2_dDelta_dTau+delta*dPSI3_dDelta2_dTau)+dDELTAbi_dTau*(2*dPSI_dDelta+delta*dPSI2_dDelta2);
        long double Line2 = 2*dDELTAbi_dDelta*(dPSI_dTau+delta*dPSI2_dDelta_dTau)+2*dDELTAbi2_dDelta_dTau*(PSI+delta*dPSI_dDelta);
        long double Line3 = dDELTAbi2_dDelta2*delta*dPSI_dTau + dDELTAbi3_dDelta2_dTau*delta*PSI;
        s[i] = ni*(Line1+Line2+Line3);
	}
	return std::accumulate(s.begin(), s.end(), 0.0);
}
long double ResidualHelmholtzNonAnalytic::dTau3(const long double &tau, const long double &delta)
{
	long double log_tau = log(tau), log_delta = log(delta);
	for (unsigned int i=0; i<N; ++i)
	{
        ResidualHelmholtzNonAnalyticElement &el = elements[i];
        long double ni = el.n, ai = el.a, bi = el.b, betai = el.beta;
        long double Ai = el.A, Bi = el.B, Ci = el.C, Di = el.D;
        long double theta=(1.0-tau)+Ai*pow(pow(delta-1.0,2),1.0/(2.0*betai));
        long double DELTA=pow(theta,2)+Bi*pow(pow(delta-1.0,2),ai);
        long double PSI=exp(-Ci*pow(delta-1.0,2)-Di*pow(tau-1.0,2));
        long double dPSI_dTau=-2.0*Di*(tau-1.0)*PSI;
        long double dDELTAbi_dTau=-2.0*theta*bi*pow(DELTA,bi-1.0);
        long double dPSI2_dTau2=(2.0*Di*pow(tau-1.0,2)-1.0)*2.0*Di*PSI;
        long double dDELTAbi2_dTau2=2.0*bi*pow(DELTA,bi-1.0)+4.0*pow(theta,2)*bi*(bi-1.0)*pow(DELTA,bi-2.0);
        long double dPSI3_dTau3=2.0*Di*PSI*(-4*Di*Di*pow(tau-1,3)+6*Di*(tau-1));
        long double dDELTAbi3_dTau3 = -12.0*theta*bi*(bi-1.0)*pow(DELTA,bi-2)-8*pow(theta,3)*bi*(bi-1)*(bi-2)*pow(DELTA,bi-3.0);

        s[i] = ni*delta*(dDELTAbi3_dTau3*PSI+(3.0*dDELTAbi2_dTau2)*dPSI_dTau+(3*dDELTAbi_dTau )*dPSI2_dTau2+pow(DELTA,bi)*dPSI3_dTau3);
	}
	return std::accumulate(s.begin(), s.end(), 0.0);
}

////
////void ResidualHelmholtzSAFTAssociating::to_json(rapidjson::Value &el, rapidjson::Document &doc)
////{
////    el.AddMember("type","phir_SAFT_associating",doc.GetAllocator());
////    el.AddMember("a",a,doc.GetAllocator());
////    el.AddMember("m",m,doc.GetAllocator());
////    el.AddMember("epsilonbar",epsilonbar,doc.GetAllocator());
////    el.AddMember("vbarn",vbarn,doc.GetAllocator());
////    el.AddMember("kappabar",kappabar,doc.GetAllocator());
////}
////double ResidualHelmholtzSAFTAssociating::Deltabar(double tau, double delta)
////{
////    return this->g(this->eta(delta))*(exp(this->epsilonbar*tau)-1)*this->kappabar;
////}   
////double ResidualHelmholtzSAFTAssociating::dDeltabar_ddelta__consttau(double tau, double delta)
////{
////    return this->dg_deta(this->eta(delta))*(exp(this->epsilonbar*tau)-1)*this->kappabar*this->vbarn;
////}
////double ResidualHelmholtzSAFTAssociating::d2Deltabar_ddelta2__consttau(double tau, double delta)
////{
////    return this->d2g_deta2(this->eta(delta))*(exp(this->epsilonbar*tau)-1)*this->kappabar*pow(this->vbarn,(int)2);
////}
////double ResidualHelmholtzSAFTAssociating::dDeltabar_dtau__constdelta(double tau, double delta)
////{
////    return this->g(this->eta(delta))*this->kappabar*exp(this->epsilonbar*tau)*this->epsilonbar;
////}
////double ResidualHelmholtzSAFTAssociating::d2Deltabar_dtau2__constdelta(double tau, double delta)
////{
////    return this->g(this->eta(delta))*this->kappabar*exp(this->epsilonbar*tau)*pow(this->epsilonbar,(int)2);
////}
////double ResidualHelmholtzSAFTAssociating::d2Deltabar_ddelta_dtau(double tau, double delta)
////{
////    return this->dg_deta(this->eta(delta))*exp(this->epsilonbar*tau)*this->epsilonbar*this->kappabar*this->vbarn;
////}
////double ResidualHelmholtzSAFTAssociating::d3Deltabar_dtau3__constdelta(double tau, double delta)
////{
////    return this->g(this->eta(delta))*this->kappabar*exp(this->epsilonbar*tau)*pow(this->epsilonbar,(int)3);
////}
////double ResidualHelmholtzSAFTAssociating::d3Deltabar_ddelta_dtau2(double tau, double delta)
////{
////    return this->dg_deta(this->eta(delta))*this->kappabar*exp(this->epsilonbar*tau)*pow(this->epsilonbar,(int)2)*this->vbarn;
////}
////double ResidualHelmholtzSAFTAssociating::d3Deltabar_ddelta2_dtau(double tau, double delta)
////{
////    return this->d2g_deta2(this->eta(delta))*exp(this->epsilonbar*tau)*this->epsilonbar*this->kappabar*pow(this->vbarn,(int)2);
////}
////double ResidualHelmholtzSAFTAssociating::d3Deltabar_ddelta3__consttau(double tau, double delta)
////{
////    return this->d3g_deta3(this->eta(delta))*(exp(this->epsilonbar*tau)-1)*this->kappabar*pow(this->vbarn,(int)3);
////}
////
////double ResidualHelmholtzSAFTAssociating::X(double delta, double Deltabar)
////{
////    return 2/(sqrt(1+4*Deltabar*delta)+1);
////}
////double ResidualHelmholtzSAFTAssociating::dX_dDeltabar__constdelta(double delta, double Deltabar)
////{
////    double X = this->X(delta,Deltabar);
////    return -delta*X*X/(2*Deltabar*delta*X+1);
////}
////double ResidualHelmholtzSAFTAssociating::dX_ddelta__constDeltabar(double delta, double Deltabar)
////{
////    double X = this->X(delta,Deltabar);
////    return -Deltabar*X*X/(2*Deltabar*delta*X+1);
////}
////double ResidualHelmholtzSAFTAssociating::dX_dtau(double tau, double delta)
////{
////    double Deltabar = this->Deltabar(tau, delta);
////    return this->dX_dDeltabar__constdelta(delta, Deltabar)*this->dDeltabar_dtau__constdelta(tau, delta);
////}
////double ResidualHelmholtzSAFTAssociating::dX_ddelta(double tau, double delta)
////{
////    double Deltabar = this->Deltabar(tau, delta);
////    return (this->dX_ddelta__constDeltabar(delta, Deltabar)
////           + this->dX_dDeltabar__constdelta(delta, Deltabar)*this->dDeltabar_ddelta__consttau(tau, delta));
////}
////double ResidualHelmholtzSAFTAssociating::d2X_dtau2(double tau, double delta)
////{
////    double Deltabar = this->Deltabar(tau, delta);
////    double X = this->X(delta, Deltabar);
////    double beta = this->dDeltabar_dtau__constdelta(tau, delta);
////    double d_dXdtau_dbeta = -delta*X*X/(2*Deltabar*delta*X+1);
////    double d_dXdtau_dDeltabar = 2*delta*delta*X*X*X/pow(2*Deltabar*delta*X+1,2)*beta;
////    double d_dXdtau_dX = -2*beta*delta*X*(Deltabar*delta*X+1)/pow(2*Deltabar*delta*X+1,2);
////    double dbeta_dtau = this->d2Deltabar_dtau2__constdelta(tau, delta);
////    double dX_dDeltabar = this->dX_dDeltabar__constdelta(delta, Deltabar);
////    return d_dXdtau_dX*dX_dDeltabar*beta+d_dXdtau_dDeltabar*beta+d_dXdtau_dbeta*dbeta_dtau;
////}
////double ResidualHelmholtzSAFTAssociating::d2X_ddeltadtau(double tau, double delta)
////{
////    double Deltabar = this->Deltabar(tau, delta);
////    double X = this->X(delta, Deltabar);
////    double alpha = this->dDeltabar_ddelta__consttau(tau, delta);
////    double beta = this->dDeltabar_dtau__constdelta(tau, delta);
////    double dalpha_dtau = this->d2Deltabar_ddelta_dtau(tau, delta);
////    double d_dXddelta_dDeltabar = X*X*(2*delta*delta*X*alpha-1)/pow(2*Deltabar*delta*X+1,2);
////    double d_dXddelta_dalpha = -delta*X*X/(2*Deltabar*delta*X+1);
////    double d_dXddelta_dX = -(Deltabar+delta*alpha)*2*(Deltabar*delta*X*X+X)/pow(2*Deltabar*delta*X+1,2);
////    double dX_dDeltabar = this->dX_dDeltabar__constdelta(delta, Deltabar);
////    return d_dXddelta_dX*dX_dDeltabar*beta+d_dXddelta_dDeltabar*beta+d_dXddelta_dalpha*dalpha_dtau;
////}
////double ResidualHelmholtzSAFTAssociating::d2X_ddelta2(double tau, double delta)
////{
////    double Deltabar = this->Deltabar(tau, delta);
////    double X = this->X(delta, Deltabar);
////    double alpha = this->dDeltabar_ddelta__consttau(tau, delta);
////    double dalpha_ddelta = this->d2Deltabar_ddelta2__consttau(tau, delta);
////    
////    double dX_ddelta_constall = X*X*(2*Deltabar*Deltabar*X-alpha)/pow(2*Deltabar*delta*X+1,2);
////    double d_dXddelta_dX = -(Deltabar+delta*alpha)*2*(Deltabar*delta*X*X+X)/pow(2*Deltabar*delta*X+1,2);
////    double d_dXddelta_dDeltabar = X*X*(2*delta*delta*X*alpha-1)/pow(2*Deltabar*delta*X+1,2);
////    double d_dXddelta_dalpha = -delta*X*X/(2*Deltabar*delta*X+1);
////    
////    double dX_dDeltabar = this->dX_dDeltabar__constdelta(delta, Deltabar);
////    double dX_ddelta = this->dX_ddelta__constDeltabar(delta, Deltabar);
////
////    double val = (dX_ddelta_constall
////            + d_dXddelta_dX*dX_ddelta
////            + d_dXddelta_dX*dX_dDeltabar*alpha
////            + d_dXddelta_dDeltabar*alpha
////            + d_dXddelta_dalpha*dalpha_ddelta);
////    return val;
////}   
////double ResidualHelmholtzSAFTAssociating::d3X_dtau3(double tau, double delta)
////{
////    double Delta = this->Deltabar(tau, delta);
////    double X = this->X(delta, Delta);
////    double dX_dDelta = this->dX_dDeltabar__constdelta(delta, Delta);
////    double Delta_t = this->dDeltabar_dtau__constdelta(tau, delta);
////    double Delta_tt = this->d2Deltabar_dtau2__constdelta(tau, delta);
////    double Delta_ttt = this->d3Deltabar_dtau3__constdelta(tau, delta);
////    double dXtt_dX = 2*X*delta*(-6*Delta*pow(Delta_t, 2)*pow(X, 2)*pow(delta, 2)*(Delta*X*delta + 1) + 3*pow(Delta_t, 2)*X*delta*(2*Delta*X*delta + 1) - Delta_tt*pow(2*Delta*X*delta + 1, 3) + X*delta*(Delta*Delta_tt + 3*pow(Delta_t, 2))*pow(2*Delta*X*delta + 1, 2))/pow(2*Delta*X*delta + 1, 4);
////    double dXtt_dDelta = 2*pow(X, 3)*pow(delta, 2)*(-6*pow(Delta_t, 2)*X*delta*(Delta*X*delta + 1) - 3*pow(Delta_t, 2)*X*delta*(2*Delta*X*delta + 1) + Delta_tt*pow(2*Delta*X*delta + 1, 2))/pow(2*Delta*X*delta + 1, 4);
////    double dXtt_dDelta_t = 4*Delta_t*pow(X, 3)*pow(delta, 2)*(3*Delta*X*delta + 2)/pow(2*Delta*X*delta + 1, 3);
////    double dXtt_dDelta_tt = -pow(X, 2)*delta/(2*Delta*X*delta + 1);
////    return dXtt_dX*dX_dDelta*Delta_t+dXtt_dDelta*Delta_t + dXtt_dDelta_t*Delta_tt + dXtt_dDelta_tt*Delta_ttt;
////}
////double ResidualHelmholtzSAFTAssociating::d3X_ddeltadtau2(double tau, double delta)
////{
////    double Delta = this->Deltabar(tau, delta);
////    double X = this->X(delta, Delta);
////    double dX_ddelta = this->dX_ddelta__constDeltabar(delta, Delta);
////    double dX_dDelta = this->dX_dDeltabar__constdelta(delta, Delta);
////    double Delta_t = this->dDeltabar_dtau__constdelta(tau, delta);
////    double Delta_d = this->dDeltabar_ddelta__consttau(tau, delta);
////    double Delta_dt = this->d2Deltabar_ddelta_dtau(tau, delta);
////    double Delta_tt = this->d2Deltabar_dtau2__constdelta(tau, delta);
////    double Delta_dtt = this->d3Deltabar_ddelta_dtau2(tau, delta);
////    double dXtt_ddelta = pow(X, 2)*(-12*Delta*pow(Delta_t, 2)*pow(X, 2)*pow(delta, 2)*(Delta*X*delta + 1) + 2*pow(Delta_t, 2)*X*delta*(-Delta*X*delta + 2)*(2*Delta*X*delta + 1) - Delta_tt*pow(2*Delta*X*delta + 1, 3) + 2*X*delta*(Delta*Delta_tt + 2*pow(Delta_t, 2))*pow(2*Delta*X*delta + 1, 2))/pow(2*Delta*X*delta + 1, 4);
////    double dXtt_dX = 2*X*delta*(-6*Delta*pow(Delta_t, 2)*pow(X, 2)*pow(delta, 2)*(Delta*X*delta + 1) + 3*pow(Delta_t, 2)*X*delta*(2*Delta*X*delta + 1) - Delta_tt*pow(2*Delta*X*delta + 1, 3) + X*delta*(Delta*Delta_tt + 3*pow(Delta_t, 2))*pow(2*Delta*X*delta + 1, 2))/pow(2*Delta*X*delta + 1, 4);
////    double dXtt_dDelta = 2*pow(X, 3)*pow(delta, 2)*(-6*pow(Delta_t, 2)*X*delta*(Delta*X*delta + 1) - 3*pow(Delta_t, 2)*X*delta*(2*Delta*X*delta + 1) + Delta_tt*pow(2*Delta*X*delta + 1, 2))/pow(2*Delta*X*delta + 1, 4);
////    double dXtt_dDelta_t = 4*Delta_t*pow(X, 3)*pow(delta, 2)*(3*Delta*X*delta + 2)/pow(2*Delta*X*delta + 1, 3);
////    double dXtt_dDelta_tt = -pow(X, 2)*delta/(2*Delta*X*delta + 1);
////    return dXtt_ddelta + dXtt_dX*dX_ddelta + dXtt_dX*dX_dDelta*Delta_d + dXtt_dDelta*Delta_d + dXtt_dDelta_t*Delta_dt + dXtt_dDelta_tt*Delta_dtt;
////}
////
////double ResidualHelmholtzSAFTAssociating::d3X_ddelta2dtau(double tau, double delta)
////{
////    double Delta = this->Deltabar(tau, delta);
////    double X = this->X(delta, Delta);
////    double dX_ddelta = this->dX_ddelta__constDeltabar(delta, Delta);
////    double dX_dDelta = this->dX_dDeltabar__constdelta(delta, Delta);
////    double Delta_t = this->dDeltabar_dtau__constdelta(tau, delta);
////    double Delta_d = this->dDeltabar_ddelta__consttau(tau, delta);
////    double Delta_dd = this->d2Deltabar_ddelta2__consttau(tau, delta);
////    double Delta_dt = this->d2Deltabar_ddelta_dtau(tau, delta);
////    double Delta_tt = this->d2Deltabar_dtau2__constdelta(tau, delta);
////    double Delta_ddt = this->d3Deltabar_ddelta2_dtau(tau, delta);
////    double dXdd_dX = 2*X*(-6*Delta*pow(X, 2)*delta*pow(Delta + Delta_d*delta, 2)*(Delta*X*delta + 1) - Delta_dd*delta*pow(2*Delta*X*delta + 1, 3) + 2*X*(2*Delta*X*delta + 1)*(-Delta*Delta_d*delta*(2*Delta_d*X*pow(delta, 2) - 1) - Delta*delta*(2*pow(Delta, 2)*X - Delta_d) + Delta*(Delta + Delta_d*delta)*(Delta*X*delta + 1) + Delta_d*delta*(Delta + Delta_d*delta)*(Delta*X*delta + 1)) + pow(2*Delta*X*delta + 1, 2)*(3*pow(Delta, 2)*X + Delta*Delta_dd*X*pow(delta, 2) + Delta*X*(Delta + Delta_d*delta) + pow(Delta_d, 2)*X*pow(delta, 2) + Delta_d*X*delta*(Delta + Delta_d*delta) + Delta_d*(2*Delta_d*X*pow(delta, 2) - 1) - Delta_d))/pow(2*Delta*X*delta + 1, 4);
////    double dXdd_ddelta = pow(X, 2)*(-24*pow(Delta, 4)*pow(X, 3)*delta - 8*pow(Delta, 3)*Delta_d*pow(X, 3)*pow(delta, 2) - 18*pow(Delta, 3)*pow(X, 2) + 8*pow(Delta, 2)*Delta_d*pow(X, 2)*delta - 4*pow(Delta, 2)*Delta_dd*pow(X, 2)*pow(delta, 2) + 10*Delta*pow(Delta_d, 2)*pow(X, 2)*pow(delta, 2) + 12*Delta*Delta_d*X - 4*Delta*Delta_dd*X*delta + 8*pow(Delta_d, 2)*X*delta - Delta_dd)/(16*pow(Delta, 4)*pow(X, 4)*pow(delta, 4) + 32*pow(Delta, 3)*pow(X, 3)*pow(delta, 3) + 24*pow(Delta, 2)*pow(X, 2)*pow(delta, 2) + 8*Delta*X*delta + 1);
////    double dXdd_dDelta = pow(X, 3)*(-8*pow(Delta, 2)*Delta_d*pow(X, 2)*pow(delta, 3) + 8*pow(Delta, 2)*Delta_dd*pow(X, 2)*pow(delta, 4) + 10*pow(Delta, 2)*X*delta - 24*Delta*pow(Delta_d, 2)*pow(X, 2)*pow(delta, 4) + 8*Delta*Delta_d*X*pow(delta, 2) + 8*Delta*Delta_dd*X*pow(delta, 3) + 8*Delta - 18*pow(Delta_d, 2)*X*pow(delta, 3) + 12*Delta_d*delta + 2*Delta_dd*pow(delta, 2))/(16*pow(Delta, 4)*pow(X, 4)*pow(delta, 4) + 32*pow(Delta, 3)*pow(X, 3)*pow(delta, 3) + 24*pow(Delta, 2)*pow(X, 2)*pow(delta, 2) + 8*Delta*X*delta + 1);
////    double dXdd_dDelta_d = 2*pow(X, 2)*(2*X*delta*(Delta + Delta_d*delta)*(Delta*X*delta + 1) + (2*Delta*X*delta + 1)*(2*Delta_d*X*pow(delta, 2) - 1))/pow(2*Delta*X*delta + 1, 3);
////    double dXdd_dDelta_dd = -pow(X, 2)*delta/(2*Delta*X*delta + 1);
////
////    return dXdd_dX*dX_dDelta*Delta_t + dXdd_dDelta*Delta_t + dXdd_dDelta_d*Delta_dt + dXdd_dDelta_dd*Delta_ddt;
////}
////
////double Xdd(double X, double delta, double Delta, double Delta_d, double Delta_dd)
////{
////    return Delta*pow(X, 2)*(2*Delta + 2*Delta_d*delta)*(Delta*pow(X, 2)*delta + X)/pow(2*Delta*X*delta + 1, 3) + Delta_d*pow(X, 2)*delta*(2*Delta + 2*Delta_d*delta)*(Delta*pow(X, 2)*delta + X)/pow(2*Delta*X*delta + 1, 3) + Delta_d*pow(X, 2)*(2*Delta_d*X*pow(delta, 2) - 1)/pow(2*Delta*X*delta + 1, 2) - Delta_dd*pow(X, 2)*delta/(2*Delta*X*delta + 1) + pow(X, 2)*(2*pow(Delta, 2)*X - Delta_d)/pow(2*Delta*X*delta + 1, 2);
////}
////
////double ResidualHelmholtzSAFTAssociating::d3X_ddelta3(double tau, double delta)
////{
////    double Delta = this->Deltabar(tau, delta);
////    double X = this->X(delta, Delta);
////    double dX_ddelta = this->dX_ddelta__constDeltabar(delta, Delta);
////    double dX_dDelta = this->dX_dDeltabar__constdelta(delta, Delta);
////    double Delta_d = this->dDeltabar_ddelta__consttau(tau, delta);
////    double Delta_dd = this->d2Deltabar_ddelta2__consttau(tau, delta);
////    double Delta_ddd = this->d3Deltabar_ddelta3__consttau(tau, delta);
////
////    double dXdd_dX = 2*X*(-6*Delta*pow(X, 2)*delta*pow(Delta + Delta_d*delta, 2)*(Delta*X*delta + 1) - Delta_dd*delta*pow(2*Delta*X*delta + 1, 3) + 2*X*(2*Delta*X*delta + 1)*(-Delta*Delta_d*delta*(2*Delta_d*X*pow(delta, 2) - 1) - Delta*delta*(2*pow(Delta, 2)*X - Delta_d) + Delta*(Delta + Delta_d*delta)*(Delta*X*delta + 1) + Delta_d*delta*(Delta + Delta_d*delta)*(Delta*X*delta + 1)) + pow(2*Delta*X*delta + 1, 2)*(3*pow(Delta, 2)*X + Delta*Delta_dd*X*pow(delta, 2) + Delta*X*(Delta + Delta_d*delta) + pow(Delta_d, 2)*X*pow(delta, 2) + Delta_d*X*delta*(Delta + Delta_d*delta) + Delta_d*(2*Delta_d*X*pow(delta, 2) - 1) - Delta_d))/pow(2*Delta*X*delta + 1, 4);
////    double dXdd_ddelta = pow(X, 2)*(-24*pow(Delta, 4)*pow(X, 3)*delta - 8*pow(Delta, 3)*Delta_d*pow(X, 3)*pow(delta, 2) - 18*pow(Delta, 3)*pow(X, 2) + 8*pow(Delta, 2)*Delta_d*pow(X, 2)*delta - 4*pow(Delta, 2)*Delta_dd*pow(X, 2)*pow(delta, 2) + 10*Delta*pow(Delta_d, 2)*pow(X, 2)*pow(delta, 2) + 12*Delta*Delta_d*X - 4*Delta*Delta_dd*X*delta + 8*pow(Delta_d, 2)*X*delta - Delta_dd)/(16*pow(Delta, 4)*pow(X, 4)*pow(delta, 4) + 32*pow(Delta, 3)*pow(X, 3)*pow(delta, 3) + 24*pow(Delta, 2)*pow(X, 2)*pow(delta, 2) + 8*Delta*X*delta + 1);
////    double dXdd_dDelta = pow(X, 3)*(-8*pow(Delta, 2)*Delta_d*pow(X, 2)*pow(delta, 3) + 8*pow(Delta, 2)*Delta_dd*pow(X, 2)*pow(delta, 4) + 10*pow(Delta, 2)*X*delta - 24*Delta*pow(Delta_d, 2)*pow(X, 2)*pow(delta, 4) + 8*Delta*Delta_d*X*pow(delta, 2) + 8*Delta*Delta_dd*X*pow(delta, 3) + 8*Delta - 18*pow(Delta_d, 2)*X*pow(delta, 3) + 12*Delta_d*delta + 2*Delta_dd*pow(delta, 2))/(16*pow(Delta, 4)*pow(X, 4)*pow(delta, 4) + 32*pow(Delta, 3)*pow(X, 3)*pow(delta, 3) + 24*pow(Delta, 2)*pow(X, 2)*pow(delta, 2) + 8*Delta*X*delta + 1);
////    double dXdd_dDelta_d = 2*pow(X, 2)*(2*X*delta*(Delta + Delta_d*delta)*(Delta*X*delta + 1) + (2*Delta*X*delta + 1)*(2*Delta_d*X*pow(delta, 2) - 1))/pow(2*Delta*X*delta + 1, 3);
////    double dXdd_dDelta_dd = -pow(X, 2)*delta/(2*Delta*X*delta + 1);
////
////    return dXdd_ddelta + dXdd_dX*(dX_ddelta + dX_dDelta*Delta_d) + dXdd_dDelta*Delta_d + dXdd_dDelta_d*Delta_dd + dXdd_dDelta_dd*Delta_ddd;
////}
////
////
////double ResidualHelmholtzSAFTAssociating::g(double eta)
////{
////    return 0.5*(2-eta)/pow(1-eta,(int)3);
////}    
////double ResidualHelmholtzSAFTAssociating::dg_deta(double eta)
////{
////    return 0.5*(5-2*eta)/pow(1-eta,(int)4);
////}
////double ResidualHelmholtzSAFTAssociating::d2g_deta2(double eta)
////{
////    return 3*(3-eta)/pow(1-eta,(int)5);
////}   
////double ResidualHelmholtzSAFTAssociating::d3g_deta3(double eta)
////{
////    return 6*(7-2*eta)/pow(1-eta,(int)6);
////}   
////double ResidualHelmholtzSAFTAssociating::eta(double delta){
////    return this->vbarn*delta;
////}
////double ResidualHelmholtzSAFTAssociating::A(double log_tau, double tau, double log_delta, double delta, int i)
////{
////    double X = this->X(delta, this->Deltabar(tau, delta));
////    return this->m*this->a*((log(X)-X/2.0+0.5));
////}
////double ResidualHelmholtzSAFTAssociating::dA_dDelta(double log_tau, double tau, double log_delta, double delta, int i)
////{
////    double X = this->X(delta, this->Deltabar(tau, delta));
////    return this->m*this->a*(1/X-0.5)*this->dX_ddelta(tau, delta);
////}
////double ResidualHelmholtzSAFTAssociating::dA_dTau(double log_tau, double tau, double log_delta, double delta, int i)
////{
////    double X = this->X(delta, this->Deltabar(tau, delta));
////    return this->m*this->a*(1/X-0.5)*this->dX_dtau(tau, delta);
////}
////double ResidualHelmholtzSAFTAssociating::d2A_dTau2(double log_tau, double tau, double log_delta, double delta, int i)
////{
////    double X = this->X(delta, this->Deltabar(tau, delta));
////    double X_tau = this->dX_dtau(tau, delta);
////    double X_tautau = this->d2X_dtau2(tau, delta);
////    return this->m*this->a*((1/X-0.5)*X_tautau-pow(X_tau/X, 2));
////}
////double ResidualHelmholtzSAFTAssociating::d2A_dDelta2(double log_tau, double tau, double log_delta, double delta, int i)
////{
////    double X = this->X(delta, this->Deltabar(tau, delta));
////    double X_delta = this->dX_ddelta(tau, delta);
////    double X_deltadelta = this->d2X_ddelta2(tau, delta);
////    return this->m*this->a*((1/X-0.5)*X_deltadelta-pow(X_delta/X,2));
////}
////double ResidualHelmholtzSAFTAssociating::d2A_dDelta_dTau(double log_tau, double tau, double log_delta, double delta, int i)
////{
////    double X = this->X(delta, this->Deltabar(tau, delta));
////    double X_delta = this->dX_ddelta(tau, delta);
////    double X_deltadelta = this->d2X_ddelta2(tau, delta);
////    double X_tau = this->dX_dtau(tau, delta);
////    double X_deltatau = this->d2X_ddeltadtau(tau, delta);
////    return this->m*this->a*((-X_tau/X/X)*X_delta+X_deltatau*(1/X-0.5));
////}
////double ResidualHelmholtzSAFTAssociating::d3A_dTau3(double log_tau, double tau, double log_delta, double delta, int i)
////{
////    double X = this->X(delta, this->Deltabar(tau, delta));
////    double X_t = this->dX_dtau(tau, delta);
////    double X_tt = this->d2X_dtau2(tau, delta);
////    double X_ttt = this->d3X_dtau3(tau, delta);
////    return this->m*this->a*((1/X-1.0/2.0)*X_ttt+(-X_t/pow(X,(int)2))*X_tt-2*(pow(X,(int)2)*(X_t*X_tt)-pow(X_t,(int)2)*(X*X_t))/pow(X,(int)4));
////}
////double ResidualHelmholtzSAFTAssociating::d3A_dDelta_dTau2(double log_tau, double tau, double log_delta, double delta, int i)
////{
////    double X = this->X(delta, this->Deltabar(tau, delta));
////    double X_t = this->dX_dtau(tau, delta);
////    double X_d = this->dX_ddelta(tau, delta);
////    double X_tt = this->d2X_dtau2(tau, delta);
////    double X_dt = this->d2X_ddeltadtau(tau, delta);
////    double X_dtt = this->d3X_ddeltadtau2(tau, delta);
////    return this->m*this->a*((1/X-1.0/2.0)*X_dtt-X_d/pow(X,(int)2)*X_tt-2*(pow(X,(int)2)*(X_t*X_dt)-pow(X_t,(int)2)*(X*X_d))/pow(X,(int)4));
////}
////double ResidualHelmholtzSAFTAssociating::d3A_dDelta2_dTau(double log_tau, double tau, double log_delta, double delta, int i)
////{
////    double X = this->X(delta, this->Deltabar(tau, delta));
////    double X_t = this->dX_dtau(tau, delta);
////    double X_d = this->dX_ddelta(tau, delta);
////    double X_dd = this->d2X_ddelta2(tau, delta);
////    double X_dt = this->d2X_ddeltadtau(tau, delta);
////    double X_ddt = this->d3X_ddelta2dtau(tau, delta);
////    return this->m*this->a*((1/X-1.0/2.0)*X_ddt-X_t/pow(X,(int)2)*X_dd-2*(pow(X,(int)2)*(X_d*X_dt)-pow(X_d,(int)2)*(X*X_t))/pow(X,(int)4));
////}
////double ResidualHelmholtzSAFTAssociating::d3A_dDelta3(double log_tau, double tau, double log_delta, double delta, int i)
////{
////    double X = this->X(delta, this->Deltabar(tau, delta));
////    double X_d = this->dX_ddelta(tau, delta);
////    double X_dd = this->d2X_ddelta2(tau, delta);
////    double X_ddd = this->d3X_ddelta3(tau, delta);
////    return this->m*this->a*((1/X-1.0/2.0)*X_ddd-X_d/pow(X,(int)2)*X_dd-2*(pow(X,(int)2)*(X_d*X_dd)-pow(X_d,(int)2)*(X*X_d))/pow(X,(int)4));
////}

}; /* namespace CoolProp */