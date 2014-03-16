#include <vector>
#include "rapidjson/rapidjson_include.h"
#include "Helmholtz.h"


std::vector<double> ResidualHelmholtzTerm::dDelta(std::vector<double> tau, std::vector<double> delta)
{
    std::vector<double> out = tau;
    for (int i = 0; i < (int)tau.size(); i++)
    {
        out[i] = dDelta(tau[i],delta[i]);
    }
    return out;
}
std::vector<double> ResidualHelmholtzTerm::dDelta2(std::vector<double> tau, std::vector<double> delta)
{
    std::vector<double> out = tau;
    for (int i = 0; i < (int)tau.size(); i++)
    {
        out[i] = dDelta2(tau[i],delta[i]);
    }
    return out;
}
std::vector<double> ResidualHelmholtzTerm::dTau2(std::vector<double> tau, std::vector<double> delta)
{
    std::vector<double> out = tau;
    for (int i = 0; i < (int)tau.size(); i++)
    {
        out[i] = dTau2(tau[i],delta[i]);
    }
    return out;
}
std::vector<double> ResidualHelmholtzTerm::dDelta_dTau(std::vector<double> tau, std::vector<double> delta)
{
    std::vector<double> out = tau;
    for (int i = 0; i < (int)tau.size(); i++)
    {
        out[i] = dDelta_dTau(tau[i],delta[i]);
    }
    return out;
}

void ResidualHelmholtzPower::to_json(rapidjson::Value &el, rapidjson::Document &doc)
{
    el.AddMember("type","alphar_power",doc.GetAllocator());
    rapidjson::Value _n(rapidjson::kArrayType), _d(rapidjson::kArrayType), _t(rapidjson::kArrayType), _l(rapidjson::kArrayType);
    for (unsigned int i = 0; i <= n.size(); ++i)
    {
        _n.PushBack(n[i],doc.GetAllocator());
        _d.PushBack(d[i],doc.GetAllocator());
        _t.PushBack(t[i],doc.GetAllocator());
        _l.PushBack(l[i],doc.GetAllocator());
    }
    el.AddMember("n",_n,doc.GetAllocator());
    el.AddMember("d",_d,doc.GetAllocator());
    el.AddMember("t",_t,doc.GetAllocator());
    el.AddMember("l",_l,doc.GetAllocator());
}

double ResidualHelmholtzTerm::eval(int key, double tau, double delta)
{
    double summer=0, log_tau = log(tau), log_delta = log(delta);
    switch (key)
    {
    case idA_dDelta:
        { for (unsigned int i = 0; i <= n.size(); ++i) { summer += n[i]*dA_dDelta(log_tau, tau, log_delta, delta, i); } break; }
    case idA_dTau:
        { for (unsigned int i = 0; i <= n.size(); ++i) { summer += n[i]*dA_dTau(log_tau, tau, log_delta, delta, i); } break; }
    case id2A_dDelta_dTau:
        { for (unsigned int i = 0; i <= n.size(); ++i) { summer += n[i]*d2A_dDelta_dTau(log_tau, tau, log_delta, delta, i); } break; }
    };
    return summer;
}
double ResidualHelmholtzPower::A(double log_tau, double tau, double log_delta, double delta, int i)
{
    if (l[i]>0)
        return exp(t[i]*log_tau+d[i]*log_delta-pow(delta,(int)l[i]));
    else
        return exp(t[i]*log_tau+d[i]*log_delta);
}
double ResidualHelmholtzPower::dA_dDelta(double log_tau, double tau, double log_delta, double delta, int i)
{
    double pow_delta_li, li, ni, di, ti;
    ni = n[i]; di = d[i]; ti = t[i]; li = l[i]; 
    if (li > 0){
        pow_delta_li = pow(delta,(int)li);
        return (di-li*pow_delta_li)*exp(ti*log_tau+(di-1)*log_delta-pow_delta_li);
    }
    else
    {
        return di*exp(ti*log_tau+(di-1)*log_delta);
    }
}
double ResidualHelmholtzPower::dA_dTau(double log_tau, double tau, double log_delta, double delta, int i)
{
    if (l[i]>0)
        return t[i]*exp((t[i]-1)*log_tau+d[i]*log_delta-pow(delta,(int)l[i]));
    else
        return t[i]*exp((t[i]-1)*log_tau+d[i]*log_delta);
}
double ResidualHelmholtzPower::d2A_dDelta2(double log_tau, double tau, double log_delta, double delta, int i)
{
    if (l[i]>0){
        double pow_delta_li = pow(delta,(int)l[i]);
        return ((d[i]-l[i]*pow_delta_li)*(d[i]-1.0-l[i]*pow_delta_li) - l[i]*l[i]*pow_delta_li)*exp(t[i]*log_tau+(d[i]-2)*log_delta-pow_delta_li);
    }
    else
        return d[i]*(d[i]-1.0)*exp(t[i]*log_tau+(d[i]-2)*log_delta);
}
double ResidualHelmholtzPower::d2A_dDelta_dTau(double log_tau, double tau, double log_delta, double delta, int i)
{	
    if (l[i]>0){
        double pow_delta_li = pow(delta,(int)l[i]);
        return t[i]*(d[i]-l[i]*pow_delta_li)*exp((t[i]-1)*log_tau+(d[i]-1)*log_delta-pow_delta_li);
    }
    else
        return d[i]*t[i]*exp((t[i]-1)*log_tau+(d[i]-1)*log_delta);
}
double ResidualHelmholtzPower::d2A_dTau2(double log_tau, double tau, double log_delta, double delta, int i)
{
    if (l[i]>0)
        return t[i]*(t[i]-1)*exp((t[i]-2)*log_tau+d[i]*log_delta-pow(delta,(int)l[i]));
    else
        return t[i]*(t[i]-1)*exp((t[i]-2)*log_tau+d[i]*log_delta);
}
double ResidualHelmholtzPower::d3A_dDelta3(double log_tau, double tau, double log_delta, double delta, int i)
{
    if (l[i]>0)	
    {
        double pow_delta_li = pow(delta,(int)l[i]);
        double bracket = (d[i]*(d[i]-1)*(d[i]-2)+pow_delta_li*(-2*l[i]+6*d[i]*l[i]-3*d[i]*d[i]*l[i]-3*d[i]*l[i]*l[i]+3*l[i]*l[i]-l[i]*l[i]*l[i])+pow_delta_li*pow_delta_li*(3*d[i]*l[i]*l[i]-3*l[i]*l[i]+3*l[i]*l[i]*l[i])-l[i]*l[i]*l[i]*pow_delta_li*pow_delta_li*pow_delta_li);
        return bracket*exp(t[i]*log_tau+(d[i]-3)*log_delta-pow_delta_li);
    }
    else
        return d[i]*(d[i]-1.0)*(d[i]-2)*exp(t[i]*log_tau+(d[i]-3)*log_delta);
}
double ResidualHelmholtzPower::d3A_dDelta_dTau2(double log_tau, double tau, double log_delta, double delta, int i)
{
    if (l[i]>0){
        double pow_delta_li = pow(delta,(int)l[i]);
        return t[i]*(t[i]-1)*(d[i]-l[i]*pow_delta_li)*exp((t[i]-2)*log_tau+(d[i]-1)*log_delta-pow_delta_li);
    }
    else
        return t[i]*(t[i]-1)*d[i]*exp((t[i]-2)*log_tau+(d[i]-1)*log_delta);
}

double ResidualHelmholtzPower::d3A_dDelta2_dTau(double log_tau, double tau, double log_delta, double delta, int i)
{
    if (l[i]>0){
        double pow_delta_li = pow(delta,(int)l[i]);
        return t[i]*(((d[i]-l[i]*pow_delta_li))*(d[i]-1-l[i]*pow_delta_li)-l[i]*l[i]*pow_delta_li)*exp((t[i]-1)*log_tau+(d[i]-2)*log_delta-pow_delta_li);
    }
    else
        return d[i]*t[i]*(d[i]-1)*exp((t[i]-1)*log_tau+(d[i]-2)*log_delta);
}

double ResidualHelmholtzPower::d3A_dTau3(double log_tau, double tau, double log_delta, double delta, int i)
{
    if (l[i]>0)
        return t[i]*(t[i]-1)*(t[i]-2)*exp((t[i]-3)*log_tau+d[i]*log_delta-pow(delta,(int)l[i]));
    else
        return t[i]*(t[i]-1)*(t[i]-2)*exp((t[i]-3)*log_tau+d[i]*log_delta);
}

void ResidualHelmholtzExponential::to_json(rapidjson::Value &el, rapidjson::Document &doc)
{
    el.AddMember("type","alphar_exponential",doc.GetAllocator());
    rapidjson::Value _n(rapidjson::kArrayType), _d(rapidjson::kArrayType), _t(rapidjson::kArrayType), _l(rapidjson::kArrayType), _g(rapidjson::kArrayType);
    for (unsigned int i=0; i<=n.size(); ++i)
    {
        _n.PushBack(n[i],doc.GetAllocator());
        _d.PushBack(d[i],doc.GetAllocator());
        _t.PushBack(t[i],doc.GetAllocator());
        _l.PushBack(l[i],doc.GetAllocator());
        _g.PushBack(g[i],doc.GetAllocator());
    }
    el.AddMember("n",_n,doc.GetAllocator());
    el.AddMember("d",_d,doc.GetAllocator());
    el.AddMember("t",_t,doc.GetAllocator());
    el.AddMember("l",_l,doc.GetAllocator());
    el.AddMember("g",_g,doc.GetAllocator());
}

// Term and its derivatives
double ResidualHelmholtzExponential::A(double log_tau, double tau, double log_delta, double delta, int i)
{
    return exp(t[i]*log_tau+d[i]*log_delta-g[i]*pow(delta,l[i]));
}
double ResidualHelmholtzExponential::dA_dTau(double log_tau, double tau, double log_delta, double delta, int i)
{
    return t[i]*exp((t[i]-1)*log_tau+d[i]*log_delta-g[i]*pow(delta,l[i]));
}
double ResidualHelmholtzExponential::d2A_dTau2(double log_tau, double tau, double log_delta, double delta, int i)
{
    return t[i]*(t[i]-1)*exp((t[i]-2)*log_tau+d[i]*log_delta-g[i]*pow(delta,l[i]));
}
double ResidualHelmholtzExponential::d3A_dTau3(double log_tau, double tau, double log_delta, double delta, int i)
{
    return t[i]*(t[i]-1)*(t[i]-2)*exp((t[i]-3)*log_tau+d[i]*log_delta-g[i]*pow(delta,l[i]));
}
double ResidualHelmholtzExponential::d3A_dDelta_dTau2(double log_tau, double tau, double log_delta, double delta, int i)
{
    double pow_delta_li = pow(delta,l[i]);
    return t[i]*(t[i]-1)*(d[i]-g[i]*l[i]*pow_delta_li)*exp((t[i]-2)*log_tau+(d[i]-1)*log_delta-g[i]*pow_delta_li);
}
double ResidualHelmholtzExponential::dA_dDelta(double log_tau, double tau, double log_delta, double delta, int i)
{
    double pow_delta_li = pow(delta,l[i]);
    return (d[i]-g[i]*l[i]*pow_delta_li)*exp(t[i]*log_tau+(d[i]-1)*log_delta-g[i]*pow_delta_li);
    
}
double ResidualHelmholtzExponential::d2A_dDelta2(double log_tau, double tau, double log_delta, double delta, int i)
{
    double pow_delta_li = pow(delta,l[i]);
    // Typo in Span, 2000, re-derived from Sympy
    double bracket = d[i]*d[i] - 2*d[i]*pow(delta,l[i])*g[i]*l[i] - d[i] + pow(delta,2*l[i])*g[i]*g[i]*l[i]*l[i] - pow(delta,l[i])*g[i]*l[i]*l[i] + pow(delta,l[i])*g[i]*l[i];
    return bracket*exp(t[i]*log_tau+(d[i]-2)*log_delta-g[i]*pow_delta_li);
}
double ResidualHelmholtzExponential::d3A_dDelta3(double log_tau, double tau, double log_delta, double delta, int i)
{
    // >> n_i, tau, t_i, d_i, delta, g_i, l_i = symbols(' n_i tau t_i d_i delta g_i l_i')
    // >> phir = n_i*tau**t_i*delta**d_i*exp(-g_i*pow(delta,l_i))
    // >> simplify(diff(diff(diff(phir,delta),delta),delta))
    double pow_delta_li = pow(delta,l[i]);
    double pow_delta_2li = pow(delta,2*l[i]);
    double pow_delta_3li = pow(delta,3*l[i]);
    double bracket = d[i]*d[i]*d[i] - 3*d[i]*d[i]*pow_delta_li*g[i]*l[i] - 3*d[i]*d[i] + 3*d[i]*pow_delta_2li*g[i]*g[i]*l[i]*l[i] - 3*d[i]*pow_delta_li*g[i]*l[i]*l[i] + 6*d[i]*pow_delta_li*g[i]*l[i] + 2*d[i] - pow_delta_3li*g[i]*g[i]*g[i]*l[i]*l[i]*l[i] + 3*pow_delta_2li*g[i]*g[i]*l[i]*l[i]*l[i] - 3*pow_delta_2li*g[i]*g[i]*l[i]*l[i] - pow_delta_li*g[i]*l[i]*l[i]*l[i] + 3*pow_delta_li*g[i]*l[i]*l[i] - 2*pow_delta_li*g[i]*l[i];
    return bracket*exp(t[i]*log_tau+(d[i]-3)*log_delta-g[i]*pow_delta_li);
}
double ResidualHelmholtzExponential::d3A_dDelta2_dTau(double log_tau, double tau, double log_delta, double delta, int i)
{
    double pow_delta_li = pow(delta,l[i]);
    // Typo in Span, 2000, re-derived from Sympy
    double bracket = d[i]*d[i] - 2*d[i]*pow(delta,l[i])*g[i]*l[i] - d[i] + pow(delta,2*l[i])*g[i]*g[i]*l[i]*l[i] - pow(delta,l[i])*g[i]*l[i]*l[i] + pow(delta,l[i])*g[i]*l[i];
    return t[i]*bracket*exp((t[i]-1)*log_tau+(d[i]-2)*log_delta-g[i]*pow_delta_li);
}
double ResidualHelmholtzExponential::d2A_dDelta_dTau(double log_tau, double tau, double log_delta, double delta, int i)
{
    double pow_delta_li = pow(delta,l[i]);
    return t[i]*(d[i]-g[i]*l[i]*pow_delta_li)*exp((t[i]-1)*log_tau+(d[i]-1)*log_delta-g[i]*pow_delta_li);
}



void ResidualHelmholtzGaussian::to_json(rapidjson::Value &el, rapidjson::Document &doc)
{
    el.AddMember("type","alphar_gaussian",doc.GetAllocator());
    rapidjson::Value _n(rapidjson::kArrayType), _d(rapidjson::kArrayType), _t(rapidjson::kArrayType), 
        _eta(rapidjson::kArrayType), _epsilon(rapidjson::kArrayType), _beta(rapidjson::kArrayType), _gamma(rapidjson::kArrayType);
    for (unsigned int i=0;i<=n.size();++i)
    {
        _n.PushBack(n[i],doc.GetAllocator());
        _d.PushBack(d[i],doc.GetAllocator());
        _t.PushBack(t[i],doc.GetAllocator());
        _eta.PushBack(eta[i],doc.GetAllocator());
        _epsilon.PushBack(epsilon[i],doc.GetAllocator());
        _beta.PushBack(beta[i],doc.GetAllocator());
        _gamma.PushBack(gamma[i],doc.GetAllocator());
    }
    el.AddMember("n",_n,doc.GetAllocator());
    el.AddMember("d",_d,doc.GetAllocator());
    el.AddMember("t",_t,doc.GetAllocator());
    el.AddMember("eta",_eta,doc.GetAllocator());
    el.AddMember("epsilon",_epsilon,doc.GetAllocator());
    el.AddMember("beta",_beta,doc.GetAllocator());
    el.AddMember("gamma",_gamma,doc.GetAllocator());
}

// Term and its derivatives
double ResidualHelmholtzGaussian::A(double log_tau, double tau, double log_delta, double delta, int i)
{
    double psi=exp(-eta[i]*pow(delta-epsilon[i],2)-beta[i]*pow(tau-gamma[i],2));
    return pow(delta,d[i])*pow(tau,t[i])*psi;
}
double ResidualHelmholtzGaussian::dA_dTau(double log_tau, double tau, double log_delta, double delta, int i)
{
    double psi=exp(-eta[i]*pow(delta-epsilon[i],2)-beta[i]*pow(tau-gamma[i],2));
    return pow(delta,d[i])*pow(tau,t[i])*psi*(t[i]/tau-2.0*beta[i]*(tau-gamma[i]));
}
double ResidualHelmholtzGaussian::d2A_dTau2(double log_tau, double tau, double log_delta, double delta, int i)
{
    double psi=exp(-eta[i]*pow(delta-epsilon[i],2)-beta[i]*pow(tau-gamma[i],2));
    return pow(delta,d[i])*pow(tau,t[i])*psi*(pow(t[i]/tau-2.0*beta[i]*(tau-gamma[i]),2)-t[i]/pow(tau,2)-2.0*beta[i]);
}
double ResidualHelmholtzGaussian::d3A_dTau3(double log_tau, double tau, double log_delta, double delta, int i)
{
    // triple derivative product rule (a*b*c)' = a'*b*c+a*b'*c+a*b*c'
    double psi=exp(-eta[i]*pow(delta-epsilon[i],2)-beta[i]*pow(tau-gamma[i],2));
    double dpsi_dTau = exp(-eta[i]*pow(delta-epsilon[i],2)-beta[i]*pow(tau-gamma[i],2))*(-2*beta[i]*(tau-gamma[i]));

    double bracket = pow(t[i]/tau-2.0*beta[i]*(tau-gamma[i]),2)-t[i]/pow(tau,2)-2.0*beta[i];
    double dbracket_dTau = 2*(t[i]/tau-2.0*beta[i]*(tau-gamma[i]))*(-t[i]/tau/tau-2*beta[i])+2*t[i]/pow(tau,3);
    return pow(delta,d[i])*(t[i]*pow(tau,t[i]-1)*psi*bracket+pow(tau,t[i])*dpsi_dTau*bracket+pow(tau,t[i])*psi*dbracket_dTau);
}
double ResidualHelmholtzGaussian::dA_dDelta(double log_tau, double tau, double log_delta, double delta, int i)
{
    double psi=exp(-eta[i]*pow(delta-epsilon[i],2)-beta[i]*pow(tau-gamma[i],2));
    return pow(delta,d[i])*pow(tau,t[i])*psi*(d[i]/delta-2.0*eta[i]*(delta-epsilon[i]));
}
double ResidualHelmholtzGaussian::d2A_dDelta2(double log_tau, double tau, double log_delta, double delta, int i)
{
    double psi=exp(-eta[i]*pow(delta-epsilon[i],2)-beta[i]*pow(tau-gamma[i],2));
    return pow(tau,t[i])*psi*(-2.0*eta[i]*pow(delta,d[i])+4.0*pow(eta[i],2)*pow(delta,d[i])*pow(delta-epsilon[i],2)-4.0*d[i]*eta[i]*pow(delta,d[i]-1)*(delta-epsilon[i])+d[i]*(d[i]-1.0)*pow(delta,d[i]-2));
}
double ResidualHelmholtzGaussian::d3A_dDelta3(double log_tau, double tau, double log_delta, double delta, int i)
{
    double psi=exp(-eta[i]*pow(delta-epsilon[i],2)-beta[i]*pow(tau-gamma[i],2));
    double bracket = (pow(d[i]-2*eta[i]*delta*(delta-epsilon[i]),3)-3*d[i]*d[i]+2*d[i]-6*d[i]*eta[i]*delta*delta+6*eta[i]*delta*(delta-epsilon[i])*(d[i]+2*eta[i]*delta*delta));
    return pow(tau,t[i])*pow(delta,d[i]-3)*psi*bracket;
}
double ResidualHelmholtzGaussian::d3A_dDelta2_dTau(double log_tau, double tau, double log_delta, double delta, int i)
{
    double psi=exp(-eta[i]*pow(delta-epsilon[i],2)-beta[i]*pow(tau-gamma[i],2));
    return pow(tau,t[i])*psi*(t[i]/tau-2.0*beta[i]*(tau-gamma[i]))*(-2.0*eta[i]*pow(delta,d[i])+4.0*pow(eta[i],2)*pow(delta,d[i])*pow(delta-epsilon[i],2)-4.0*d[i]*eta[i]*pow(delta,d[i]-1)*(delta-epsilon[i])+d[i]*(d[i]-1.0)*pow(delta,d[i]-2));
}
double ResidualHelmholtzGaussian::d2A_dDelta_dTau(double log_tau, double tau, double log_delta, double delta, int i)
{
    double psi=exp(-eta[i]*pow(delta-epsilon[i],2)-beta[i]*pow(tau-gamma[i],2));
    return pow(delta,d[i])*pow(tau,t[i])*psi*(d[i]/delta-2.0*eta[i]*(delta-epsilon[i]))*(t[i]/tau-2.0*beta[i]*(tau-gamma[i]));
}
double ResidualHelmholtzGaussian::d3A_dDelta_dTau2(double log_tau, double tau, double log_delta, double delta, int i)
{
    double psi=exp(-eta[i]*pow(delta-epsilon[i],2)-beta[i]*pow(tau-gamma[i],2));
    return pow(delta,d[i])*pow(tau,t[i])*psi*(d[i]/delta-2.0*eta[i]*(delta-epsilon[i]))*(pow(t[i]-2.0*beta[i]*tau*(tau-gamma[i]),2)-t[i]-2*beta[i]*tau*tau)/tau/tau;
}

void ResidualHelmholtzLemmon2005::to_json(rapidjson::Value &el, rapidjson::Document &doc)
{
    el.AddMember("type","alphar_Lemmon2005",doc.GetAllocator());
    rapidjson::Value _n(rapidjson::kArrayType), _d(rapidjson::kArrayType), _t(rapidjson::kArrayType), _l(rapidjson::kArrayType), _m(rapidjson::kArrayType);
    for (unsigned int i=0;i<=n.size();i++)
    {
        _n.PushBack(n[i],doc.GetAllocator());
        _d.PushBack(d[i],doc.GetAllocator());
        _t.PushBack(t[i],doc.GetAllocator());
        _l.PushBack(l[i],doc.GetAllocator());
        _m.PushBack(m[i],doc.GetAllocator());
    }
    el.AddMember("n",_n,doc.GetAllocator());
    el.AddMember("d",_d,doc.GetAllocator());
    el.AddMember("t",_t,doc.GetAllocator());
    el.AddMember("l",_l,doc.GetAllocator());
    el.AddMember("m",_m,doc.GetAllocator());
}

// Term and its derivatives
double ResidualHelmholtzLemmon2005::A(double log_tau, double tau, double log_delta, double delta, int i)
{
    if (l[i] != 0 && m[i] != 0)
        return exp(t[i]*log_tau+d[i]*log_delta-pow(delta,l[i])-pow(tau,m[i]));
    else if (l[i] != 0 && m[i] == 0)
        return exp(t[i]*log_tau+d[i]*log_delta-pow(delta,l[i]));
    else
        return exp(t[i]*log_tau+d[i]*log_delta);
}
double ResidualHelmholtzLemmon2005::dA_dTau(double log_tau, double tau, double log_delta, double delta, int i)
{
    if (l[i] != 0 && m[i] != 0){
        double pow_tau_mi = pow(tau,m[i]);
        return (t[i]-m[i]*pow_tau_mi)*exp((t[i]-1)*log_tau+d[i]*log_delta-pow(delta,l[i])-pow_tau_mi);
    }
    else if (l[i] != 0 && m[i] == 0)
        return t[i]*exp((t[i]-1)*log_tau+d[i]*log_delta-pow(delta,l[i]));
    else
        return t[i]*exp((t[i]-1)*log_tau+d[i]*log_delta);
}
double ResidualHelmholtzLemmon2005::d2A_dTau2(double log_tau, double tau, double log_delta, double delta, int i)
{
    if (l[i] != 0 && m[i] != 0){
        double pow_tau_mi = pow(tau,m[i]);
        double bracket = (t[i]-m[i]*pow_tau_mi)*(t[i]-1-m[i]*pow_tau_mi)-m[i]*m[i]*pow_tau_mi;
        return bracket*exp((t[i]-2)*log_tau+d[i]*log_delta-pow(delta,l[i])-pow_tau_mi);
    }
    else if (l[i] != 0 && m[i] == 0)
        return t[i]*(t[i]-1)*exp((t[i]-2)*log_tau+d[i]*log_delta-pow(delta,l[i]));
    else
        return t[i]*(t[i]-1)*exp((t[i]-2)*log_tau+d[i]*log_delta);
}

/*!

\f[
\frac{{{\partial ^2}{\alpha ^r}}}{{\partial {\tau ^2}}} = {N_k}{\delta ^{{d_k}}}{\tau ^{{t_k} - 2}}\exp \left( { - {\delta ^{{l_k}}}} \right)\exp \left( { - {\tau ^{{m_k}}}} \right)\left[ {\left( {{t_k} - {m_k}{\tau ^{{m_k}}}} \right)\left( {{t_k} - 1 - {m_k}{\tau ^{{m_k}}}} \right) - m_k^2{\tau ^{{m_k}}}} \right]\\
\f]
\f[
\frac{{{\partial ^2}{\alpha ^r}}}{{\partial {\tau ^2}}} = {N_k}{\delta ^{{d_k}}}\exp \left( { - {\delta ^{{l_k}}}} \right){\tau ^{{t_k} - 2}}\exp \left( { - {\tau ^{{m_k}}}} \right)\left[ {\left( {{t_k} - {m_k}{\tau ^{{m_k}}}} \right)\left( {{t_k} - 1 - {m_k}{\tau ^{{m_k}}}} \right) - m_k^2{\tau ^{{m_k}}}} \right]\\
\f]
Group all the terms that don't depend on \$ \tau \$
\f[
\frac{{{\partial ^2}{\alpha ^r}}}{{\partial {\tau ^2}}} = A{\tau ^{{t_k} - 2}}\exp \left( { - {\tau ^{{m_k}}}} \right)\left[ {\left( {{t_k} - {m_k}{\tau ^{{m_k}}}} \right)\left( {{t_k} - 1 - {m_k}{\tau ^{{m_k}}}} \right) - m_k^2{\tau ^{{m_k}}}} \right]\\
\f]
\f[
\frac{1}{A}\frac{{{\partial ^3}{\alpha ^r}}}{{\partial {\tau ^3}}} = {\tau ^{{t_k} - 2}}\exp \left( { - {\tau ^{{m_k}}}} \right)\frac{\partial }{{\partial \tau }}\left[ {\left( {{t_k} - {m_k}{\tau ^{{m_k}}}} \right)\left( {{t_k} - 1 - {m_k}{\tau ^{{m_k}}}} \right) - m_k^2{\tau ^{{m_k}}}} \right] + \frac{\partial }{{\partial \tau }}\left[ {{\tau ^{{t_k} - 2}}\exp \left( { - {\tau ^{{m_k}}}} \right)} \right]\left[ {\left( {{t_k} - {m_k}{\tau ^{{m_k}}}} \right)\left( {{t_k} - 1 - {m_k}{\tau ^{{m_k}}}} \right) - m_k^2{\tau ^{{m_k}}}} \right]\\
\f]
\f[
\frac{\partial }{{\partial \tau }}\left[ {{\tau ^{{t_k} - 2}}\exp \left( { - {\tau ^{{m_k}}}} \right)} \right] = ({t_k} - 2){\tau ^{{t_k} - 3}}\exp \left( { - {\tau ^{{m_k}}}} \right) + {\tau ^{{t_k} - 2}}\exp \left( { - {\tau ^{{m_k}}}} \right)( - {m_k}{\tau ^{{m_k} - 1}}) = \exp \left( { - {\tau ^{{m_k}}}} \right)\left( {({t_k} - 2){\tau ^{{t_k} - 3}} - {\tau ^{{t_k} - 2}}{m_k}{\tau ^{{m_k} - 1}}} \right)\\
\f]
\f[
\frac{\partial }{{\partial \tau }}\left[ {\left( {{t_k} - {m_k}{\tau ^{{m_k}}}} \right)\left( {{t_k} - 1 - {m_k}{\tau ^{{m_k}}}} \right) - m_k^2{\tau ^{{m_k}}}} \right] = \left( {{t_k} - {m_k}{\tau ^{{m_k}}}} \right)\left( { - m_k^2{\tau ^{{m_k} - 1}}} \right) + \left( { - m_k^2{\tau ^{{m_k} - 1}}} \right)\left( {{t_k} - 1 - {m_k}{\tau ^{{m_k}}}} \right) - m_k^3{\tau ^{{m_k} - 1}} =  - m_k^2{\tau ^{{m_k} - 1}}\left[ {{t_k} - {m_k}{\tau ^{{m_k}}} + {t_k} - 1 - {m_k}{\tau ^{{m_k}}} + {m_k}} \right] =  - m_k^2{\tau ^{{m_k} - 1}}\left[ {2{t_k} - 2{m_k}{\tau ^{{m_k}}} - 1 + {m_k}} \right]\\
\f]
\f[
\frac{1}{A}\frac{{{\partial ^3}{\alpha ^r}}}{{\partial {\tau ^3}}} = {\tau ^{{t_k} - 2}}\exp \left( { - {\tau ^{{m_k}}}} \right)\left( { - m_k^2{\tau ^{{m_k} - 1}}\left[ {2{t_k} - 2{m_k}{\tau ^{{m_k}}} - 1 + {m_k}} \right]} \right) + \exp \left( { - {\tau ^{{m_k}}}} \right)\left( {({t_k} - 2){\tau ^{{t_k} - 3}} - {\tau ^{{t_k} - 2}}{m_k}{\tau ^{{m_k} - 1}}} \right)\left[ {\left( {{t_k} - {m_k}{\tau ^{{m_k}}}} \right)\left( {{t_k} - 1 - {m_k}{\tau ^{{m_k}}}} \right) - m_k^2{\tau ^{{m_k}}}} \right]\\
\f]
\f[
\frac{1}{A}\frac{{{\partial ^3}{\alpha ^r}}}{{\partial {\tau ^3}}} = \exp \left( { - {\tau ^{{m_k}}}} \right)\left[ { - {\tau ^{{t_k} - 2}}m_k^2{\tau ^{{m_k} - 1}}\left[ {2{t_k} - 2{m_k}{\tau ^{{m_k}}} - 1 + {m_k}} \right] + \left( {({t_k} - 2){\tau ^{{t_k} - 3}} - {\tau ^{{t_k} - 2}}{m_k}{\tau ^{{m_k} - 1}}} \right)\left[ {\left( {{t_k} - {m_k}{\tau ^{{m_k}}}} \right)\left( {{t_k} - 1 - {m_k}{\tau ^{{m_k}}}} \right) - m_k^2{\tau ^{{m_k}}}} \right]} \right]
\f]
*/
double ResidualHelmholtzLemmon2005::d3A_dTau3(double log_tau, double tau, double log_delta, double delta, int i)
{
    if (l[i] != 0 && m[i] != 0){
        double pow_delta_li = pow(delta,l[i]);
        double pow_tau_mi = pow(tau,m[i]);
        double bracket = -pow(tau,t[i]+m[i]-3)*m[i]*m[i]*(2*t[i]-2*m[i]*pow_tau_mi-1-m[i])+((t[i]-2)*pow(tau,t[i]-3)-pow(tau,t[i]-2)*m[i]*pow(tau,m[i]-1))*((t[i]-m[i]*pow_tau_mi)*(t[i]-1-m[i]*pow_tau_mi)-m[i]*m[i]*pow_tau_mi);
        return t[i]*(t[i]-1)*(t[i]-2)*exp((t[i]-3)*log_tau+d[i]*log_delta-pow_delta_li-pow_tau_mi);
    }
    else if (l[i] != 0 && m[i] == 0){
        return t[i]*(t[i]-1)*(t[i]-2)*exp((t[i]-3)*log_tau+d[i]*log_delta-pow(delta,l[i]));
    }
    else
        return t[i]*(t[i]-1)*(t[i]-2)*exp((t[i]-3)*log_tau+d[i]*log_delta);
}
double ResidualHelmholtzLemmon2005::d3A_dDelta_dTau2(double log_tau, double tau, double log_delta, double delta, int i)
{
    if (l[i] != 0 && m[i] != 0){
        double pow_delta_li = pow(delta,l[i]);
        double pow_tau_mi = pow(tau,m[i]);
        // delta derivative of second tau derivative
        double bracket = ((t[i]-m[i]*pow_tau_mi)*(t[i]-1-m[i]*pow_tau_mi)-m[i]*m[i]*pow_tau_mi)*(d[i]-l[i]*pow_delta_li);
        return bracket*exp((t[i]-2)*log_tau+(d[i]-1)*log_delta-pow_delta_li-pow_tau_mi);
    }
    else if (l[i] != 0 && m[i] == 0){
        double pow_delta_li = pow(delta,l[i]);
        return t[i]*(t[i]-1)*(d[i]-l[i]*pow_delta_li)*exp((t[i]-2)*log_tau+(d[i]-1)*log_delta-pow_delta_li);
    }
    else
        return t[i]*(t[i]-1)*d[i]*exp((t[i]-2)*log_tau+(d[i]-1)*log_delta);
}
double ResidualHelmholtzLemmon2005::dA_dDelta(double log_tau, double tau, double log_delta, double delta, int i)
{
    if (l[i] != 0 && m[i] != 0){
        double pow_delta_li = pow(delta,l[i]);
        double pow_tau_mi = pow(tau,m[i]);
        return (d[i]-l[i]*pow_delta_li)*exp(t[i]*log_tau+(d[i]-1)*log_delta-pow_delta_li-pow_tau_mi);
    }
    else if (l[i]>0 && m[i] == 0){
        double pow_delta_li = pow(delta,l[i]);
        return (d[i]-l[i]*pow_delta_li)*exp(t[i]*log_tau+(d[i]-1)*log_delta-pow_delta_li);
    }
    else
        return d[i]*exp(t[i]*log_tau+(d[i]-1)*log_delta);
}
double ResidualHelmholtzLemmon2005::d2A_dDelta2(double log_tau, double tau, double log_delta, double delta, int i)
{
    if (l[i] != 0 && m[i] != 0){	
        double pow_delta_li = pow(delta,l[i]);
        double pow_tau_mi = pow(tau,m[i]);
        return ((d[i]-l[i]*pow_delta_li)*(d[i]-1.0-l[i]*pow_delta_li) - l[i]*l[i]*pow_delta_li)*exp(t[i]*log_tau+(d[i]-2)*log_delta-pow_delta_li-pow_tau_mi);
    }
        
    else if (l[i] != 0 && m[i] == 0){
        double pow_delta_li = pow(delta,l[i]);
        return ((d[i]-l[i]*pow_delta_li)*(d[i]-1.0-l[i]*pow_delta_li) - l[i]*l[i]*pow_delta_li)*exp(t[i]*log_tau+(d[i]-2)*log_delta-pow_delta_li);
    }
    else
        return d[i]*(d[i]-1.0)*exp(t[i]*log_tau+(d[i]-2)*log_delta);
}
double ResidualHelmholtzLemmon2005::d3A_dDelta3(double log_tau, double tau, double log_delta, double delta, int i)
{
    if (l[i] != 0 && m[i] != 0){
        double pow_delta_li = pow(delta,l[i]);
        double pow_tau_mi = pow(tau,m[i]);
        double bracket = (d[i]*(d[i]-1)*(d[i]-2)+pow_delta_li*(-2*l[i]+6*d[i]*l[i]-3*d[i]*d[i]*l[i]-3*d[i]*l[i]*l[i]+3*l[i]*l[i]-l[i]*l[i]*l[i])+pow_delta_li*pow_delta_li*(3*d[i]*l[i]*l[i]-3*l[i]*l[i]+3*l[i]*l[i]*l[i])-l[i]*l[i]*l[i]*pow_delta_li*pow_delta_li*pow_delta_li);
        return bracket*exp(t[i]*log_tau+(d[i]-3)*log_delta-pow_delta_li-pow_tau_mi);
    }
    else if (l[i] != 0 && m[i] == 0)
    {
        double pow_delta_li = pow(delta,l[i]);
        double bracket = (d[i]*(d[i]-1)*(d[i]-2)+pow_delta_li*(-2*l[i]+6*d[i]*l[i]-3*d[i]*d[i]*l[i]-3*d[i]*l[i]*l[i]+3*l[i]*l[i]-l[i]*l[i]*l[i])+pow_delta_li*pow_delta_li*(3*d[i]*l[i]*l[i]-3*l[i]*l[i]+3*l[i]*l[i]*l[i])-l[i]*l[i]*l[i]*pow_delta_li*pow_delta_li*pow_delta_li);
        return bracket*exp(t[i]*log_tau+(d[i]-3)*log_delta-pow_delta_li);
    }
    else
        return d[i]*(d[i]-1.0)*(d[i]-2)*exp(t[i]*log_tau+(d[i]-3)*log_delta);
}

double ResidualHelmholtzLemmon2005::d3A_dDelta2_dTau(double log_tau, double tau, double log_delta, double delta, int i)
{
    if (l[i] != 0 && m[i] != 0){
        double pow_delta_li = pow(delta,l[i]);
        double pow_tau_mi = pow(tau,m[i]);
        double bracket = (t[i]-m[i]*pow_tau_mi)*(((d[i]-l[i]*pow_delta_li))*(d[i]-1-l[i]*pow_delta_li)-l[i]*l[i]*pow_delta_li);
        return bracket*exp((t[i]-1)*log_tau+(d[i]-2)*log_delta-pow_delta_li-pow_tau_mi);
    }
    else if (l[i] != 0 && m[i] == 0){
        double pow_delta_li = pow(delta,l[i]);
        double bracket = t[i]*(((d[i]-l[i]*pow_delta_li))*(d[i]-1-l[i]*pow_delta_li)-l[i]*l[i]*pow_delta_li);
        return bracket*exp((t[i]-1)*log_tau+(d[i]-2)*log_delta-pow_delta_li);
    }
    else
        return d[i]*t[i]*(d[i]-1)*exp((t[i]-1)*log_tau+(d[i]-2)*log_delta);
}
double ResidualHelmholtzLemmon2005::d2A_dDelta_dTau(double log_tau, double tau, double log_delta, double delta, int i)
{
    if (l[i] != 0 && m[i] != 0){
        double pow_delta_li = pow(delta,l[i]);
        double pow_tau_mi = pow(tau,m[i]);
        return (d[i]-l[i]*pow_delta_li)*(t[i]-m[i]*pow_tau_mi)*exp((t[i]-1)*log_tau+(d[i]-1)*log_delta-pow_delta_li-pow_tau_mi);
    }
    else if (l[i] != 0 && m[i] == 0){
        double pow_delta_li = pow(delta,l[i]);
        return t[i]*(d[i]-l[i]*pow_delta_li)*exp((t[i]-1)*log_tau+(d[i]-1)*log_delta-pow_delta_li);
    }
    else
        return d[i]*t[i]*exp((t[i]-1)*log_tau+(d[i]-1)*log_delta);
}




void ResidualHelmholtzGERG2008Gaussian::to_json(rapidjson::Value &el, rapidjson::Document &doc)
{
    el.AddMember("type","alphar_GERG2008_gaussian",doc.GetAllocator());
    rapidjson::Value _n(rapidjson::kArrayType), _d(rapidjson::kArrayType), _t(rapidjson::kArrayType), 
        _eta(rapidjson::kArrayType), _epsilon(rapidjson::kArrayType), _beta(rapidjson::kArrayType), _gamma(rapidjson::kArrayType);
    for (unsigned int i=0; i<=n.size(); ++i)
    {
        _n.PushBack(n[i],doc.GetAllocator());
        _d.PushBack(d[i],doc.GetAllocator());
        _t.PushBack(t[i],doc.GetAllocator());
        _eta.PushBack(eta[i],doc.GetAllocator());
        _epsilon.PushBack(epsilon[i],doc.GetAllocator());
        _beta.PushBack(beta[i],doc.GetAllocator());
        _gamma.PushBack(gamma[i],doc.GetAllocator());
    }
    el.AddMember("n",_n,doc.GetAllocator());
    el.AddMember("d",_d,doc.GetAllocator());
    el.AddMember("t",_t,doc.GetAllocator());
    el.AddMember("eta",_eta,doc.GetAllocator());
    el.AddMember("epsilon",_epsilon,doc.GetAllocator());
    el.AddMember("beta",_beta,doc.GetAllocator());
    el.AddMember("gamma",_gamma,doc.GetAllocator());
}

double ResidualHelmholtzGERG2008Gaussian::A(double log_tau, double tau, double log_delta, double delta, int i)
{
    return pow(tau,t[i])*pow(delta,d[i])*psi(tau,delta,i);
}
double ResidualHelmholtzGERG2008Gaussian::dA_dDelta(double log_tau, double tau, double log_delta, double delta, int i)
{
    return pow(tau,t[i])*pow(delta,d[i])*psi(tau,delta,i)*(d[i]/delta-2*eta[i]*(delta-epsilon[i])-beta[i]);
}
double ResidualHelmholtzGERG2008Gaussian::d2A_dDelta2(double log_tau, double tau, double log_delta, double delta, int i)
{	
    return pow(tau,t[i])*pow(delta,d[i])*psi(tau,delta,i)*(pow(d[i]/delta-2*eta[i]*(delta-epsilon[i])-beta[i],2)-d[i]/delta/delta-2*eta[i]);
}
double ResidualHelmholtzGERG2008Gaussian::dA_dTau(double log_tau, double tau, double log_delta, double delta, int i)
{
    return t[i]*pow(tau,t[i]-1)*pow(delta,d[i])*psi(tau,delta,i);
}
double ResidualHelmholtzGERG2008Gaussian::d2A_dTau2(double log_tau, double tau, double log_delta, double delta, int i)
{
    return t[i]*(t[i]-1)*pow(tau,t[i]-2)*pow(delta,d[i])*psi(tau,delta,i);
}
double ResidualHelmholtzGERG2008Gaussian::d2A_dDelta_dTau(double log_tau, double tau, double log_delta, double delta, int i)
{
    return t[i]*pow(tau,t[i]-1)*pow(delta,d[i])*psi(tau,delta,i)*(d[i]/delta-2*eta[i]*(delta-epsilon[i])-beta[i]);
}

void ResidualHelmholtzNonAnalytic::to_json(rapidjson::Value &el, rapidjson::Document &doc)
{
    el.AddMember("type","alphar_critical",doc.GetAllocator());

    rapidjson::Value _n(rapidjson::kArrayType), _a(rapidjson::kArrayType), _b(rapidjson::kArrayType), 
        _beta(rapidjson::kArrayType), __A(rapidjson::kArrayType), _B(rapidjson::kArrayType), _C(rapidjson::kArrayType), _D(rapidjson::kArrayType);
    for (unsigned int i=0; i<=n.size(); ++i)
    {
        _n.PushBack(n[i],doc.GetAllocator());
        _a.PushBack(a[i],doc.GetAllocator());
        _b.PushBack(b[i],doc.GetAllocator());
        _beta.PushBack(beta[i],doc.GetAllocator());
        __A.PushBack(_A[i],doc.GetAllocator());
        _B.PushBack(B[i],doc.GetAllocator());
        _C.PushBack(C[i],doc.GetAllocator());
        _D.PushBack(D[i],doc.GetAllocator());
    }
    el.AddMember("n",_n,doc.GetAllocator());
    el.AddMember("a",_a,doc.GetAllocator());
    el.AddMember("b",_b,doc.GetAllocator());
    el.AddMember("beta",_beta,doc.GetAllocator());
    el.AddMember("A",__A,doc.GetAllocator());
    el.AddMember("B",_B,doc.GetAllocator());
    el.AddMember("C",_C,doc.GetAllocator());
    el.AddMember("D",_D,doc.GetAllocator());
}

double ResidualHelmholtzNonAnalytic::A(double log_tau, double tau, double log_delta, double delta, int i)
{
    double theta=(1.0-tau)+_A[i]*pow(pow(delta-1.0,2),1/(2*beta[i]));
    double DELTA=pow(theta,2)+B[i]*pow(pow(delta-1.0,2),a[i]);
    double PSI=exp(-C[i]*pow(delta-1.0,2)-D[i]*pow(tau-1.0,2));
    return pow(DELTA,b[i])*delta*PSI;
}

double ResidualHelmholtzNonAnalytic::dA_dDelta(double log_tau, double tau, double log_delta, double delta, int i)
{
    double dDELTAbi_dDelta;
    double theta=(1.0-tau)+_A[i]*pow(pow(delta-1.0,2),1.0/(2.0*beta[i]));
    double DELTA=pow(theta,2)+B[i]*pow(pow(delta-1.0,2),a[i]);
    double PSI=exp(-C[i]*pow(delta-1.0,2)-D[i]*pow(tau-1.0,2));
    double dPSI_dDelta=-2.0*C[i]*(delta-1.0)*PSI;
    double dDELTA_dDelta=(delta-1.0)*(_A[i]*theta*2.0/beta[i]*pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-1.0)+2.0*B[i]*a[i]*pow(pow(delta-1.0,2),a[i]-1.0));
        
    // At critical point, DELTA is 0, and 1/0^n is undefined
    if (fabs(DELTA) < 10*DBL_EPSILON)
    {
        dDELTAbi_dDelta = 0;
    }
    else{
        dDELTAbi_dDelta=b[i]*pow(DELTA,b[i]-1.0)*dDELTA_dDelta;
    }
    return (pow(DELTA,b[i])*(PSI+delta*dPSI_dDelta)+dDELTAbi_dDelta*delta*PSI);
}
double ResidualHelmholtzNonAnalytic::d3A_dDelta_dTau2(double log_tau, double tau, double log_delta, double delta, int i)
{
    double theta=(1.0-tau)+_A[i]*pow(pow(delta-1.0,2),1.0/(2.0*beta[i]));
    double DELTA=pow(theta,2)+B[i]*pow(pow(delta-1.0,2),a[i]);
    double PSI=exp(-C[i]*pow(delta-1.0,2)-D[i]*pow(tau-1.0,2));
        
    double dPSI_dDelta=-2.0*C[i]*(delta-1.0)*PSI;
    double dDELTA_dDelta=(delta-1.0)*(_A[i]*theta*2.0/beta[i]*pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-1.0)+2.0*B[i]*a[i]*pow(pow(delta-1.0,2),a[i]-1.0));
    double dDELTAbi_dDelta=b[i]*pow(DELTA,b[i]-1.0)*dDELTA_dDelta;
        
    double dPSI2_dDelta2=(2.0*C[i]*pow(delta-1.0,2)-1.0)*2.0*C[i]*PSI;
    double dDELTA2_dDelta2=1.0/(delta-1.0)*dDELTA_dDelta+pow(delta-1.0,2)*(4.0*B[i]*a[i]*(a[i]-1.0)*pow(pow(delta-1.0,2),a[i]-2.0)+2.0*pow(_A[i]/beta[i],2)*pow(pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-1.0),2)+_A[i]*theta*4.0/beta[i]*(1.0/(2.0*beta[i])-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-2.0));
    double dDELTAbi2_dDelta2=b[i]*(pow(DELTA,b[i]-1.0)*dDELTA2_dDelta2+(b[i]-1.0)*pow(DELTA,b[i]-2.0)*pow(dDELTA_dDelta,2));
        
    double dPSI_dTau=-2.0*D[i]*(tau-1.0)*PSI;
    double dDELTAbi_dTau=-2.0*theta*b[i]*pow(DELTA,b[i]-1.0);

    double dPSI2_dTau2=(2.0*D[i]*pow(tau-1.0,2)-1.0)*2.0*D[i]*PSI;
    double dDELTAbi2_dTau2=2.0*b[i]*pow(DELTA,b[i]-1.0)+4.0*pow(theta,2)*b[i]*(b[i]-1.0)*pow(DELTA,b[i]-2.0);
    double dPSI2_dDelta_dTau=4.0*C[i]*D[i]*(delta-1.0)*(tau-1.0)*PSI;
    double dDELTAbi2_dDelta_dTau=-_A[i]*b[i]*2.0/beta[i]*pow(DELTA,b[i]-1.0)*(delta-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-1.0)-2.0*theta*b[i]*(b[i]-1.0)*pow(DELTA,b[i]-2.0)*dDELTA_dDelta;

    double dPSI3_dDelta_dTau2 = 2*D[i]*(2*D[i]*pow(tau-1,2)-1)*dPSI_dDelta;
    double dtheta_dDelta = _A[i]/(2*beta[i])*pow(pow(delta-1,2),1/(2*beta[i])-1)*2*(delta-1);
    double dDELTAbi3_dDelta_dTau2 = 2*b[i]*(b[i]-1)*pow(DELTA,b[i]-2)*dDELTA_dDelta+4*pow(theta,2)*b[i]*(b[i]-1)*(b[i]-2)*pow(DELTA,b[i]-3)*dDELTA_dDelta+8*theta*b[i]*(b[i]-1)*pow(DELTA,b[i]-2)*dtheta_dDelta;
        
    return delta*(dDELTAbi2_dTau2*dPSI_dDelta+dDELTAbi3_dDelta_dTau2*PSI+2*dDELTAbi_dTau*dPSI2_dDelta_dTau+2.0*dDELTAbi2_dDelta_dTau*dPSI_dTau+pow(DELTA,b[i])*dPSI3_dDelta_dTau2+dDELTAbi_dDelta*dPSI2_dTau2)+n[i]*(dDELTAbi2_dTau2*PSI+2.0*dDELTAbi_dTau*dPSI_dTau+pow(DELTA,b[i])*dPSI2_dTau2);
}

double ResidualHelmholtzNonAnalytic::d2A_dDelta2(double log_tau, double tau, double log_delta, double delta, int i)
{
    double theta=(1.0-tau)+_A[i]*pow(pow(delta-1.0,2),1.0/(2.0*beta[i]));
    double DELTA=pow(theta,2)+B[i]*pow(pow(delta-1.0,2),a[i]);
    double PSI=exp(-C[i]*pow(delta-1.0,2)-D[i]*pow(tau-1.0,2));
        
    double dPSI_dDelta=-2.0*C[i]*(delta-1.0)*PSI;
    double dDELTA_dDelta=(delta-1.0)*(_A[i]*theta*2.0/beta[i]*pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-1.0)+2.0*B[i]*a[i]*pow(pow(delta-1.0,2),a[i]-1.0));
    double dDELTAbi_dDelta=b[i]*pow(DELTA,b[i]-1.0)*dDELTA_dDelta;
        
    double dPSI2_dDelta2=(2.0*C[i]*pow(delta-1.0,2)-1.0)*2.0*C[i]*PSI;
    double dDELTA2_dDelta2=1.0/(delta-1.0)*dDELTA_dDelta+pow(delta-1.0,2)*(4.0*B[i]*a[i]*(a[i]-1.0)*pow(pow(delta-1.0,2),a[i]-2.0)+2.0*pow(_A[i]/beta[i],2)*pow(pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-1.0),2)+_A[i]*theta*4.0/beta[i]*(1.0/(2.0*beta[i])-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-2.0));
    double dDELTAbi2_dDelta2=b[i]*(pow(DELTA,b[i]-1.0)*dDELTA2_dDelta2+(b[i]-1.0)*pow(DELTA,b[i]-2.0)*pow(dDELTA_dDelta,2));
        
    return (pow(DELTA,b[i])*(2.0*dPSI_dDelta+delta*dPSI2_dDelta2)+2.0*dDELTAbi_dDelta*(PSI+delta*dPSI_dDelta)+dDELTAbi2_dDelta2*delta*PSI);
}
double ResidualHelmholtzNonAnalytic::d3A_dDelta3(double log_tau, double tau, double log_delta, double delta, int i)
{
    double theta=(1.0-tau)+_A[i]*pow(pow(delta-1.0,2),1.0/(2.0*beta[i]));
    double DELTA=pow(theta,2)+B[i]*pow(pow(delta-1.0,2),a[i]);
    double PSI=exp(-C[i]*pow(delta-1.0,2)-D[i]*pow(tau-1.0,2));
        
    double dPSI_dDelta=-2.0*C[i]*(delta-1.0)*PSI;
    double dDELTA_dDelta=(delta-1.0)*(_A[i]*theta*2.0/beta[i]*pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-1.0)+2.0*B[i]*a[i]*pow(pow(delta-1.0,2),a[i]-1.0));
    double dDELTAbi_dDelta=b[i]*pow(DELTA,b[i]-1.0)*dDELTA_dDelta;
        
    double dPSI2_dDelta2=(2.0*C[i]*pow(delta-1.0,2)-1.0)*2.0*C[i]*PSI;
    double dDELTA2_dDelta2=1.0/(delta-1.0)*dDELTA_dDelta+pow(delta-1.0,2)*(4.0*B[i]*a[i]*(a[i]-1.0)*pow(pow(delta-1.0,2),a[i]-2.0)+2.0*pow(_A[i]/beta[i],2)*pow(pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-1.0),2)+_A[i]*theta*4.0/beta[i]*(1.0/(2.0*beta[i])-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-2.0));
    double dDELTAbi2_dDelta2=b[i]*(pow(DELTA,b[i]-1.0)*dDELTA2_dDelta2+(b[i]-1.0)*pow(DELTA,b[i]-2.0)*pow(dDELTA_dDelta,2));

    double dPSI3_dDelta3=2.0*C[i]*PSI*(-4*C[i]*C[i]*pow(delta-1.0,3)+6*C[i]*(delta-1));
    double dtheta_dDelta = _A[i]/(2*beta[i])*pow(pow(delta-1,2),1/(2*beta[i])-1)*2*(delta-1);

    double PI = 4*B[i]*a[i]*(a[i]-1)*pow(pow(delta-1,2),a[i]-2)+2*pow(_A[i]/beta[i],2)*pow(pow(delta-1,2),1/beta[i]-2)+4*_A[i]*theta/beta[i]*(1/(2*beta[i])-1)*pow(pow(delta-1,2),1/(2*beta[i])-2);
    double dPI_dDelta = -8*B[i]*a[i]*(a[i]-1)*(a[i]-2)*pow(pow(delta-1,2),a[i]-5.0/2.0)-8*pow(_A[i]/beta[i],2)*(1/(2*beta[i])-1)*pow(pow(delta-1,2),1/beta[i]-5.0/2.0)-(8*_A[i]*theta)/beta[i]*(1/(2*beta[i])-1)*(1/(2*beta[i])-2)*pow(pow(delta-1,2),1/(2*beta[i])-5.0/2.0)+4*_A[i]/beta[i]*(1/(2*beta[i])-1)*pow(pow(delta-1,2),1/(2*beta[i])-2)*dtheta_dDelta;
    double dDELTA3_dDelta3 = 1/(delta-1)*dDELTA2_dDelta2-1/pow(delta-1,2)*dDELTA_dDelta+pow(delta-1,2)*dPI_dDelta+2*(delta-1)*PI;
    double dDELTAbi3_dDelta3 = b[i]*(pow(DELTA,b[i]-1)*dDELTA3_dDelta3+dDELTA2_dDelta2*(b[i]-1)*pow(DELTA,b[i]-2)*dDELTA_dDelta+(b[i]-1)*(pow(DELTA,b[i]-2)*2*dDELTA_dDelta*dDELTA2_dDelta2+pow(dDELTA_dDelta,2)*(b[i]-2)*pow(DELTA,b[i]-3)*dDELTA_dDelta));
        
    return (pow(DELTA,b[i])*(3.0*dPSI2_dDelta2+delta*dPSI3_dDelta3)+3.0*dDELTAbi_dDelta*(2*dPSI_dDelta+delta*dPSI2_dDelta2)+3*dDELTAbi2_dDelta2*(PSI+delta*dPSI_dDelta)+dDELTAbi3_dDelta3*PSI*delta);
}

double ResidualHelmholtzNonAnalytic::d3A_dDelta2_dTau(double log_tau, double tau, double log_delta, double delta, int i)
{
    double theta=(1.0-tau)+_A[i]*pow(pow(delta-1.0,2),1.0/(2.0*beta[i]));
    double DELTA=pow(theta,2)+B[i]*pow(pow(delta-1.0,2),a[i]);
    double PSI=exp(-C[i]*pow(delta-1.0,2)-D[i]*pow(tau-1.0,2));
        
    double dPSI_dDelta=-2.0*C[i]*(delta-1.0)*PSI;
    double dDELTA_dDelta=(delta-1.0)*(_A[i]*theta*2.0/beta[i]*pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-1.0)+2.0*B[i]*a[i]*pow(pow(delta-1.0,2),a[i]-1.0));
    double dDELTAbi_dDelta=b[i]*pow(DELTA,b[i]-1.0)*dDELTA_dDelta;
        
    double dPSI2_dDelta2=(2.0*C[i]*pow(delta-1.0,2)-1.0)*2.0*C[i]*PSI;
    double dDELTA2_dDelta2=1.0/(delta-1.0)*dDELTA_dDelta+pow(delta-1.0,2)*(4.0*B[i]*a[i]*(a[i]-1.0)*pow(pow(delta-1.0,2),a[i]-2.0)+2.0*pow(_A[i]/beta[i],2)*pow(pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-1.0),2)+_A[i]*theta*4.0/beta[i]*(1.0/(2.0*beta[i])-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-2.0));
    double dDELTAbi2_dDelta2=b[i]*(pow(DELTA,b[i]-1.0)*dDELTA2_dDelta2+(b[i]-1.0)*pow(DELTA,b[i]-2.0)*pow(dDELTA_dDelta,2));
        
    double dPSI_dTau=-2.0*D[i]*(tau-1.0)*PSI;
    double dDELTAbi_dTau=-2.0*theta*b[i]*pow(DELTA,b[i]-1.0);
        
    double dPSI2_dDelta_dTau=4.0*C[i]*D[i]*(delta-1.0)*(tau-1.0)*PSI;
    double dDELTAbi2_dDelta_dTau=-_A[i]*b[i]*2.0/beta[i]*pow(DELTA,b[i]-1.0)*(delta-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-1.0)-2.0*theta*b[i]*(b[i]-1.0)*pow(DELTA,b[i]-2.0)*dDELTA_dDelta;
    
    //Following Terms added for this derivative
    double dPSI3_dDelta2_dTau = (2.0*C[i]*pow(delta-1.0,2)-1.0)*2.0*C[i]*dPSI_dTau;
    double dDELTA_dTau = -2*theta;
    double dDELTA2_dDelta_dTau = 2.0*_A[i]/(beta[i])*pow(pow(delta-1,2),1.0/(2.0*beta[i])-0.5);
    double dDELTA3_dDelta2_dTau = 2.0*_A[i]*(beta[i]-1)/(beta[i]*beta[i])*pow(pow(delta-1,2),1/(2*beta[i])-1.0);
        
    double dDELTAbim1_dTau = (b[i]-1)*pow(DELTA,b[i]-2)*dDELTA_dTau;
    double dDELTAbim2_dTau = (b[i]-2)*pow(DELTA,b[i]-3)*dDELTA_dTau;
    double Line11 = dDELTAbim1_dTau*dDELTA2_dDelta2 + pow(DELTA,b[i]-1)*dDELTA3_dDelta2_dTau;
    double Line21 = (b[i]-1)*(dDELTAbim2_dTau*pow(dDELTA_dDelta,2)+pow(DELTA,b[i]-2)*2*dDELTA_dDelta*dDELTA2_dDelta_dTau);
    double dDELTAbi3_dDelta2_dTau = b[i]*(Line11+Line21);
        
    double Line1 = pow(DELTA,b[i])*(2*dPSI2_dDelta_dTau+delta*dPSI3_dDelta2_dTau)+dDELTAbi_dTau*(2*dPSI_dDelta+delta*dPSI2_dDelta2);
    double Line2 = 2*dDELTAbi_dDelta*(dPSI_dTau+delta*dPSI2_dDelta_dTau)+2*dDELTAbi2_dDelta_dTau*(PSI+delta*dPSI_dDelta);
    double Line3 = dDELTAbi2_dDelta2*delta*dPSI_dTau + dDELTAbi3_dDelta2_dTau*delta*PSI;
    return (Line1+Line2+Line3);
}

double ResidualHelmholtzNonAnalytic::d2A_dDelta_dTau(double log_tau, double tau, double log_delta, double delta, int i)
{
    double theta=(1.0-tau)+_A[i]*pow(pow(delta-1.0,2),1.0/(2.0*beta[i]));
    double DELTA=pow(theta,2)+B[i]*pow(pow(delta-1.0,2),a[i]);
    double PSI=exp(-C[i]*pow(delta-1.0,2)-D[i]*pow(tau-1.0,2));
        
    double dPSI_dDelta=-2.0*C[i]*(delta-1.0)*PSI;
    double dDELTA_dDelta=(delta-1.0)*(_A[i]*theta*2.0/beta[i]*pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-1.0)+2.0*B[i]*a[i]*pow(pow(delta-1.0,2),a[i]-1.0));
    double dDELTAbi_dDelta=b[i]*pow(DELTA,b[i]-1.0)*dDELTA_dDelta;
    double dPSI_dTau=-2.0*D[i]*(tau-1.0)*PSI;
    double dDELTAbi_dTau=-2.0*theta*b[i]*pow(DELTA,b[i]-1.0);
        
    double dPSI2_dDelta_dTau=4.0*C[i]*D[i]*(delta-1.0)*(tau-1.0)*PSI;
    double dDELTAbi2_dDelta_dTau=-_A[i]*b[i]*2.0/beta[i]*pow(DELTA,b[i]-1.0)*(delta-1.0)*pow(pow(delta-1.0,2),1.0/(2.0*beta[i])-1.0)-2.0*theta*b[i]*(b[i]-1.0)*pow(DELTA,b[i]-2.0)*dDELTA_dDelta;
        
    return (pow(DELTA,b[i])*(dPSI_dTau+delta*dPSI2_dDelta_dTau)+delta*dDELTAbi_dDelta*dPSI_dTau+ dDELTAbi_dTau*(PSI+delta*dPSI_dDelta)+dDELTAbi2_dDelta_dTau*delta*PSI);
}
    
double ResidualHelmholtzNonAnalytic::dA_dTau(double log_tau, double tau, double log_delta, double delta, int i)
{
    double dDELTAbi_dTau;
    double theta=(1.0-tau)+_A[i]*pow(pow(delta-1.0,2),1.0/(2.0*beta[i]));
    double DELTA=pow(theta,2)+B[i]*pow(pow(delta-1.0,2),a[i]);
    double PSI=exp(-C[i]*pow(delta-1.0,2)-D[i]*pow(tau-1.0,2));
    double dPSI_dTau=-2.0*D[i]*(tau-1.0)*PSI;
    if (fabs(DELTA)<10*DBL_EPSILON)
        dDELTAbi_dTau = 0;
    else
        dDELTAbi_dTau = -2.0*theta*b[i]*pow(DELTA,b[i]-1.0);
    return delta*(dDELTAbi_dTau*PSI+pow(DELTA,b[i])*dPSI_dTau);
}

double ResidualHelmholtzNonAnalytic::d2A_dTau2(double log_tau, double tau, double log_delta, double delta, int i)
{
    double theta=(1.0-tau)+_A[i]*pow(pow(delta-1.0,2),1/(2*beta[i]));
    double DELTA=pow(theta,2)+B[i]*pow(pow(delta-1.0,2),a[i]);
    double PSI=exp(-C[i]*pow(delta-1.0,2)-D[i]*pow(tau-1.0,2));
    double dPSI_dTau=-2.0*D[i]*(tau-1.0)*PSI;
    double dDELTAbi_dTau=-2.0*theta*b[i]*pow(DELTA,b[i]-1.0);
    double dPSI2_dTau2=(2.0*D[i]*pow(tau-1.0,2)-1.0)*2.0*D[i]*PSI;
    double dDELTAbi2_dTau2=2.0*b[i]*pow(DELTA,b[i]-1.0)+4.0*pow(theta,2)*b[i]*(b[i]-1.0)*pow(DELTA,b[i]-2.0);
    return delta*(dDELTAbi2_dTau2*PSI+2.0*dDELTAbi_dTau*dPSI_dTau+pow(DELTA,b[i])*dPSI2_dTau2);
}
double ResidualHelmholtzNonAnalytic::d3A_dTau3(double log_tau, double tau, double log_delta, double delta, int i)
{
    double theta=(1.0-tau)+_A[i]*pow(pow(delta-1.0,2),1/(2*beta[i]));
    double DELTA=pow(theta,2)+B[i]*pow(pow(delta-1.0,2),a[i]);
    double PSI=exp(-C[i]*pow(delta-1.0,2)-D[i]*pow(tau-1.0,2));
    double dPSI_dTau=-2.0*D[i]*(tau-1.0)*PSI;
    double dDELTAbi_dTau=-2.0*theta*b[i]*pow(DELTA,b[i]-1.0);
    double dPSI2_dTau2=(2.0*D[i]*pow(tau-1.0,2)-1.0)*2.0*D[i]*PSI;
    double dDELTAbi2_dTau2=2.0*b[i]*pow(DELTA,b[i]-1.0)+4.0*pow(theta,2)*b[i]*(b[i]-1.0)*pow(DELTA,b[i]-2.0);
        
    double dPSI3_dTau3=2.0*D[i]*PSI*(-4*D[i]*D[i]*pow(tau-1,3)+6*D[i]*(tau-1));
    double dDELTAbi3_dTau3 = -12.0*theta*b[i]*(b[i]-1.0)*pow(DELTA,b[i]-2)-8*pow(theta,3)*b[i]*(b[i]-1)*(b[i]-2)*pow(DELTA,b[i]-3.0);
    return delta*(dDELTAbi3_dTau3*PSI+(3.0*dDELTAbi2_dTau2)*dPSI_dTau+(3*dDELTAbi_dTau )*dPSI2_dTau2+pow(DELTA,b[i])*dPSI3_dTau3);
}


void ResidualHelmholtzSAFTAssociating::to_json(rapidjson::Value &el, rapidjson::Document &doc)
{
    el.AddMember("type","phir_SAFT_associating",doc.GetAllocator());
    el.AddMember("a",a,doc.GetAllocator());
    el.AddMember("m",m,doc.GetAllocator());
    el.AddMember("epsilonbar",epsilonbar,doc.GetAllocator());
    el.AddMember("vbarn",vbarn,doc.GetAllocator());
    el.AddMember("kappabar",kappabar,doc.GetAllocator());
}
double ResidualHelmholtzSAFTAssociating::Deltabar(double tau, double delta)
{
	return this->g(this->eta(delta))*(exp(this->epsilonbar*tau)-1)*this->kappabar;
}   
double ResidualHelmholtzSAFTAssociating::dDeltabar_ddelta__consttau(double tau, double delta)
{
    return this->dg_deta(this->eta(delta))*(exp(this->epsilonbar*tau)-1)*this->kappabar*this->vbarn;
}
double ResidualHelmholtzSAFTAssociating::d2Deltabar_ddelta2__consttau(double tau, double delta)
{
    return this->d2g_deta2(this->eta(delta))*(exp(this->epsilonbar*tau)-1)*this->kappabar*pow(this->vbarn,(int)2);
}
double ResidualHelmholtzSAFTAssociating::dDeltabar_dtau__constdelta(double tau, double delta)
{
    return this->g(this->eta(delta))*this->kappabar*exp(this->epsilonbar*tau)*this->epsilonbar;
}
double ResidualHelmholtzSAFTAssociating::d2Deltabar_dtau2__constdelta(double tau, double delta)
{
    return this->g(this->eta(delta))*this->kappabar*exp(this->epsilonbar*tau)*pow(this->epsilonbar,(int)2);
}
double ResidualHelmholtzSAFTAssociating::d2Deltabar_ddelta_dtau(double tau, double delta)
{
    return this->dg_deta(this->eta(delta))*exp(this->epsilonbar*tau)*this->epsilonbar*this->kappabar*this->vbarn;
}
double ResidualHelmholtzSAFTAssociating::d3Deltabar_dtau3__constdelta(double tau, double delta)
{
    return this->g(this->eta(delta))*this->kappabar*exp(this->epsilonbar*tau)*pow(this->epsilonbar,(int)3);
}
double ResidualHelmholtzSAFTAssociating::d3Deltabar_ddelta_dtau2(double tau, double delta)
{
    return this->dg_deta(this->eta(delta))*this->kappabar*exp(this->epsilonbar*tau)*pow(this->epsilonbar,(int)2)*this->vbarn;
}
double ResidualHelmholtzSAFTAssociating::d3Deltabar_ddelta2_dtau(double tau, double delta)
{
    return this->d2g_deta2(this->eta(delta))*exp(this->epsilonbar*tau)*this->epsilonbar*this->kappabar*pow(this->vbarn,(int)2);
}
double ResidualHelmholtzSAFTAssociating::d3Deltabar_ddelta3__consttau(double tau, double delta)
{
    return this->d3g_deta3(this->eta(delta))*(exp(this->epsilonbar*tau)-1)*this->kappabar*pow(this->vbarn,(int)3);
}

double ResidualHelmholtzSAFTAssociating::X(double delta, double Deltabar)
{
	return 2/(sqrt(1+4*Deltabar*delta)+1);
}
double ResidualHelmholtzSAFTAssociating::dX_dDeltabar__constdelta(double delta, double Deltabar)
{
    double X = this->X(delta,Deltabar);
    return -delta*X*X/(2*Deltabar*delta*X+1);
}
double ResidualHelmholtzSAFTAssociating::dX_ddelta__constDeltabar(double delta, double Deltabar)
{
    double X = this->X(delta,Deltabar);
    return -Deltabar*X*X/(2*Deltabar*delta*X+1);
}
double ResidualHelmholtzSAFTAssociating::dX_dtau(double tau, double delta)
{
    double Deltabar = this->Deltabar(tau, delta);
    return this->dX_dDeltabar__constdelta(delta, Deltabar)*this->dDeltabar_dtau__constdelta(tau, delta);
}
double ResidualHelmholtzSAFTAssociating::dX_ddelta(double tau, double delta)
{
    double Deltabar = this->Deltabar(tau, delta);
    return (this->dX_ddelta__constDeltabar(delta, Deltabar)
           + this->dX_dDeltabar__constdelta(delta, Deltabar)*this->dDeltabar_ddelta__consttau(tau, delta));
}
double ResidualHelmholtzSAFTAssociating::d2X_dtau2(double tau, double delta)
{
    double Deltabar = this->Deltabar(tau, delta);
    double X = this->X(delta, Deltabar);
    double beta = this->dDeltabar_dtau__constdelta(tau, delta);
    double d_dXdtau_dbeta = -delta*X*X/(2*Deltabar*delta*X+1);
    double d_dXdtau_dDeltabar = 2*delta*delta*X*X*X/pow(2*Deltabar*delta*X+1,2)*beta;
    double d_dXdtau_dX = -2*beta*delta*X*(Deltabar*delta*X+1)/pow(2*Deltabar*delta*X+1,2);
    double dbeta_dtau = this->d2Deltabar_dtau2__constdelta(tau, delta);
    double dX_dDeltabar = this->dX_dDeltabar__constdelta(delta, Deltabar);
    return d_dXdtau_dX*dX_dDeltabar*beta+d_dXdtau_dDeltabar*beta+d_dXdtau_dbeta*dbeta_dtau;
}
double ResidualHelmholtzSAFTAssociating::d2X_ddeltadtau(double tau, double delta)
{
    double Deltabar = this->Deltabar(tau, delta);
    double X = this->X(delta, Deltabar);
    double alpha = this->dDeltabar_ddelta__consttau(tau, delta);
    double beta = this->dDeltabar_dtau__constdelta(tau, delta);
    double dalpha_dtau = this->d2Deltabar_ddelta_dtau(tau, delta);
    double d_dXddelta_dDeltabar = X*X*(2*delta*delta*X*alpha-1)/pow(2*Deltabar*delta*X+1,2);
    double d_dXddelta_dalpha = -delta*X*X/(2*Deltabar*delta*X+1);
    double d_dXddelta_dX = -(Deltabar+delta*alpha)*2*(Deltabar*delta*X*X+X)/pow(2*Deltabar*delta*X+1,2);
    double dX_dDeltabar = this->dX_dDeltabar__constdelta(delta, Deltabar);
    return d_dXddelta_dX*dX_dDeltabar*beta+d_dXddelta_dDeltabar*beta+d_dXddelta_dalpha*dalpha_dtau;
}
double ResidualHelmholtzSAFTAssociating::d2X_ddelta2(double tau, double delta)
{
    double Deltabar = this->Deltabar(tau, delta);
    double X = this->X(delta, Deltabar);
    double alpha = this->dDeltabar_ddelta__consttau(tau, delta);
    double dalpha_ddelta = this->d2Deltabar_ddelta2__consttau(tau, delta);
    
    double dX_ddelta_constall = X*X*(2*Deltabar*Deltabar*X-alpha)/pow(2*Deltabar*delta*X+1,2);
    double d_dXddelta_dX = -(Deltabar+delta*alpha)*2*(Deltabar*delta*X*X+X)/pow(2*Deltabar*delta*X+1,2);
    double d_dXddelta_dDeltabar = X*X*(2*delta*delta*X*alpha-1)/pow(2*Deltabar*delta*X+1,2);
    double d_dXddelta_dalpha = -delta*X*X/(2*Deltabar*delta*X+1);
    
    double dX_dDeltabar = this->dX_dDeltabar__constdelta(delta, Deltabar);
    double dX_ddelta = this->dX_ddelta__constDeltabar(delta, Deltabar);

    double val = (dX_ddelta_constall
            + d_dXddelta_dX*dX_ddelta
            + d_dXddelta_dX*dX_dDeltabar*alpha
            + d_dXddelta_dDeltabar*alpha
            + d_dXddelta_dalpha*dalpha_ddelta);
    return val;
}   
double ResidualHelmholtzSAFTAssociating::d3X_dtau3(double tau, double delta)
{
    double Delta = this->Deltabar(tau, delta);
    double X = this->X(delta, Delta);
    double dX_dDelta = this->dX_dDeltabar__constdelta(delta, Delta);
    double Delta_t = this->dDeltabar_dtau__constdelta(tau, delta);
    double Delta_tt = this->d2Deltabar_dtau2__constdelta(tau, delta);
    double Delta_ttt = this->d3Deltabar_dtau3__constdelta(tau, delta);
    double dXtt_dX = 2*X*delta*(-6*Delta*pow(Delta_t, 2)*pow(X, 2)*pow(delta, 2)*(Delta*X*delta + 1) + 3*pow(Delta_t, 2)*X*delta*(2*Delta*X*delta + 1) - Delta_tt*pow(2*Delta*X*delta + 1, 3) + X*delta*(Delta*Delta_tt + 3*pow(Delta_t, 2))*pow(2*Delta*X*delta + 1, 2))/pow(2*Delta*X*delta + 1, 4);
    double dXtt_dDelta = 2*pow(X, 3)*pow(delta, 2)*(-6*pow(Delta_t, 2)*X*delta*(Delta*X*delta + 1) - 3*pow(Delta_t, 2)*X*delta*(2*Delta*X*delta + 1) + Delta_tt*pow(2*Delta*X*delta + 1, 2))/pow(2*Delta*X*delta + 1, 4);
    double dXtt_dDelta_t = 4*Delta_t*pow(X, 3)*pow(delta, 2)*(3*Delta*X*delta + 2)/pow(2*Delta*X*delta + 1, 3);
    double dXtt_dDelta_tt = -pow(X, 2)*delta/(2*Delta*X*delta + 1);
    return dXtt_dX*dX_dDelta*Delta_t+dXtt_dDelta*Delta_t + dXtt_dDelta_t*Delta_tt + dXtt_dDelta_tt*Delta_ttt;
}
double ResidualHelmholtzSAFTAssociating::d3X_ddeltadtau2(double tau, double delta)
{
    double Delta = this->Deltabar(tau, delta);
    double X = this->X(delta, Delta);
    double dX_ddelta = this->dX_ddelta__constDeltabar(delta, Delta);
    double dX_dDelta = this->dX_dDeltabar__constdelta(delta, Delta);
    double Delta_t = this->dDeltabar_dtau__constdelta(tau, delta);
    double Delta_d = this->dDeltabar_ddelta__consttau(tau, delta);
    double Delta_dt = this->d2Deltabar_ddelta_dtau(tau, delta);
    double Delta_tt = this->d2Deltabar_dtau2__constdelta(tau, delta);
    double Delta_dtt = this->d3Deltabar_ddelta_dtau2(tau, delta);
    double dXtt_ddelta = pow(X, 2)*(-12*Delta*pow(Delta_t, 2)*pow(X, 2)*pow(delta, 2)*(Delta*X*delta + 1) + 2*pow(Delta_t, 2)*X*delta*(-Delta*X*delta + 2)*(2*Delta*X*delta + 1) - Delta_tt*pow(2*Delta*X*delta + 1, 3) + 2*X*delta*(Delta*Delta_tt + 2*pow(Delta_t, 2))*pow(2*Delta*X*delta + 1, 2))/pow(2*Delta*X*delta + 1, 4);
    double dXtt_dX = 2*X*delta*(-6*Delta*pow(Delta_t, 2)*pow(X, 2)*pow(delta, 2)*(Delta*X*delta + 1) + 3*pow(Delta_t, 2)*X*delta*(2*Delta*X*delta + 1) - Delta_tt*pow(2*Delta*X*delta + 1, 3) + X*delta*(Delta*Delta_tt + 3*pow(Delta_t, 2))*pow(2*Delta*X*delta + 1, 2))/pow(2*Delta*X*delta + 1, 4);
    double dXtt_dDelta = 2*pow(X, 3)*pow(delta, 2)*(-6*pow(Delta_t, 2)*X*delta*(Delta*X*delta + 1) - 3*pow(Delta_t, 2)*X*delta*(2*Delta*X*delta + 1) + Delta_tt*pow(2*Delta*X*delta + 1, 2))/pow(2*Delta*X*delta + 1, 4);
    double dXtt_dDelta_t = 4*Delta_t*pow(X, 3)*pow(delta, 2)*(3*Delta*X*delta + 2)/pow(2*Delta*X*delta + 1, 3);
    double dXtt_dDelta_tt = -pow(X, 2)*delta/(2*Delta*X*delta + 1);
    return dXtt_ddelta + dXtt_dX*dX_ddelta + dXtt_dX*dX_dDelta*Delta_d + dXtt_dDelta*Delta_d + dXtt_dDelta_t*Delta_dt + dXtt_dDelta_tt*Delta_dtt;
}

double ResidualHelmholtzSAFTAssociating::d3X_ddelta2dtau(double tau, double delta)
{
    double Delta = this->Deltabar(tau, delta);
    double X = this->X(delta, Delta);
    double dX_ddelta = this->dX_ddelta__constDeltabar(delta, Delta);
    double dX_dDelta = this->dX_dDeltabar__constdelta(delta, Delta);
    double Delta_t = this->dDeltabar_dtau__constdelta(tau, delta);
    double Delta_d = this->dDeltabar_ddelta__consttau(tau, delta);
    double Delta_dd = this->d2Deltabar_ddelta2__consttau(tau, delta);
    double Delta_dt = this->d2Deltabar_ddelta_dtau(tau, delta);
    double Delta_tt = this->d2Deltabar_dtau2__constdelta(tau, delta);
    double Delta_ddt = this->d3Deltabar_ddelta2_dtau(tau, delta);
    double dXdd_dX = 2*X*(-6*Delta*pow(X, 2)*delta*pow(Delta + Delta_d*delta, 2)*(Delta*X*delta + 1) - Delta_dd*delta*pow(2*Delta*X*delta + 1, 3) + 2*X*(2*Delta*X*delta + 1)*(-Delta*Delta_d*delta*(2*Delta_d*X*pow(delta, 2) - 1) - Delta*delta*(2*pow(Delta, 2)*X - Delta_d) + Delta*(Delta + Delta_d*delta)*(Delta*X*delta + 1) + Delta_d*delta*(Delta + Delta_d*delta)*(Delta*X*delta + 1)) + pow(2*Delta*X*delta + 1, 2)*(3*pow(Delta, 2)*X + Delta*Delta_dd*X*pow(delta, 2) + Delta*X*(Delta + Delta_d*delta) + pow(Delta_d, 2)*X*pow(delta, 2) + Delta_d*X*delta*(Delta + Delta_d*delta) + Delta_d*(2*Delta_d*X*pow(delta, 2) - 1) - Delta_d))/pow(2*Delta*X*delta + 1, 4);
    double dXdd_ddelta = pow(X, 2)*(-24*pow(Delta, 4)*pow(X, 3)*delta - 8*pow(Delta, 3)*Delta_d*pow(X, 3)*pow(delta, 2) - 18*pow(Delta, 3)*pow(X, 2) + 8*pow(Delta, 2)*Delta_d*pow(X, 2)*delta - 4*pow(Delta, 2)*Delta_dd*pow(X, 2)*pow(delta, 2) + 10*Delta*pow(Delta_d, 2)*pow(X, 2)*pow(delta, 2) + 12*Delta*Delta_d*X - 4*Delta*Delta_dd*X*delta + 8*pow(Delta_d, 2)*X*delta - Delta_dd)/(16*pow(Delta, 4)*pow(X, 4)*pow(delta, 4) + 32*pow(Delta, 3)*pow(X, 3)*pow(delta, 3) + 24*pow(Delta, 2)*pow(X, 2)*pow(delta, 2) + 8*Delta*X*delta + 1);
    double dXdd_dDelta = pow(X, 3)*(-8*pow(Delta, 2)*Delta_d*pow(X, 2)*pow(delta, 3) + 8*pow(Delta, 2)*Delta_dd*pow(X, 2)*pow(delta, 4) + 10*pow(Delta, 2)*X*delta - 24*Delta*pow(Delta_d, 2)*pow(X, 2)*pow(delta, 4) + 8*Delta*Delta_d*X*pow(delta, 2) + 8*Delta*Delta_dd*X*pow(delta, 3) + 8*Delta - 18*pow(Delta_d, 2)*X*pow(delta, 3) + 12*Delta_d*delta + 2*Delta_dd*pow(delta, 2))/(16*pow(Delta, 4)*pow(X, 4)*pow(delta, 4) + 32*pow(Delta, 3)*pow(X, 3)*pow(delta, 3) + 24*pow(Delta, 2)*pow(X, 2)*pow(delta, 2) + 8*Delta*X*delta + 1);
    double dXdd_dDelta_d = 2*pow(X, 2)*(2*X*delta*(Delta + Delta_d*delta)*(Delta*X*delta + 1) + (2*Delta*X*delta + 1)*(2*Delta_d*X*pow(delta, 2) - 1))/pow(2*Delta*X*delta + 1, 3);
    double dXdd_dDelta_dd = -pow(X, 2)*delta/(2*Delta*X*delta + 1);

    return dXdd_dX*dX_dDelta*Delta_t + dXdd_dDelta*Delta_t + dXdd_dDelta_d*Delta_dt + dXdd_dDelta_dd*Delta_ddt;
}

double Xdd(double X, double delta, double Delta, double Delta_d, double Delta_dd)
{
    return Delta*pow(X, 2)*(2*Delta + 2*Delta_d*delta)*(Delta*pow(X, 2)*delta + X)/pow(2*Delta*X*delta + 1, 3) + Delta_d*pow(X, 2)*delta*(2*Delta + 2*Delta_d*delta)*(Delta*pow(X, 2)*delta + X)/pow(2*Delta*X*delta + 1, 3) + Delta_d*pow(X, 2)*(2*Delta_d*X*pow(delta, 2) - 1)/pow(2*Delta*X*delta + 1, 2) - Delta_dd*pow(X, 2)*delta/(2*Delta*X*delta + 1) + pow(X, 2)*(2*pow(Delta, 2)*X - Delta_d)/pow(2*Delta*X*delta + 1, 2);
}

double ResidualHelmholtzSAFTAssociating::d3X_ddelta3(double tau, double delta)
{
    double Delta = this->Deltabar(tau, delta);
    double X = this->X(delta, Delta);
    double dX_ddelta = this->dX_ddelta__constDeltabar(delta, Delta);
    double dX_dDelta = this->dX_dDeltabar__constdelta(delta, Delta);
    double Delta_d = this->dDeltabar_ddelta__consttau(tau, delta);
    double Delta_dd = this->d2Deltabar_ddelta2__consttau(tau, delta);
    double Delta_ddd = this->d3Deltabar_ddelta3__consttau(tau, delta);

    double dXdd_dX = 2*X*(-6*Delta*pow(X, 2)*delta*pow(Delta + Delta_d*delta, 2)*(Delta*X*delta + 1) - Delta_dd*delta*pow(2*Delta*X*delta + 1, 3) + 2*X*(2*Delta*X*delta + 1)*(-Delta*Delta_d*delta*(2*Delta_d*X*pow(delta, 2) - 1) - Delta*delta*(2*pow(Delta, 2)*X - Delta_d) + Delta*(Delta + Delta_d*delta)*(Delta*X*delta + 1) + Delta_d*delta*(Delta + Delta_d*delta)*(Delta*X*delta + 1)) + pow(2*Delta*X*delta + 1, 2)*(3*pow(Delta, 2)*X + Delta*Delta_dd*X*pow(delta, 2) + Delta*X*(Delta + Delta_d*delta) + pow(Delta_d, 2)*X*pow(delta, 2) + Delta_d*X*delta*(Delta + Delta_d*delta) + Delta_d*(2*Delta_d*X*pow(delta, 2) - 1) - Delta_d))/pow(2*Delta*X*delta + 1, 4);
    double dXdd_ddelta = pow(X, 2)*(-24*pow(Delta, 4)*pow(X, 3)*delta - 8*pow(Delta, 3)*Delta_d*pow(X, 3)*pow(delta, 2) - 18*pow(Delta, 3)*pow(X, 2) + 8*pow(Delta, 2)*Delta_d*pow(X, 2)*delta - 4*pow(Delta, 2)*Delta_dd*pow(X, 2)*pow(delta, 2) + 10*Delta*pow(Delta_d, 2)*pow(X, 2)*pow(delta, 2) + 12*Delta*Delta_d*X - 4*Delta*Delta_dd*X*delta + 8*pow(Delta_d, 2)*X*delta - Delta_dd)/(16*pow(Delta, 4)*pow(X, 4)*pow(delta, 4) + 32*pow(Delta, 3)*pow(X, 3)*pow(delta, 3) + 24*pow(Delta, 2)*pow(X, 2)*pow(delta, 2) + 8*Delta*X*delta + 1);
    double dXdd_dDelta = pow(X, 3)*(-8*pow(Delta, 2)*Delta_d*pow(X, 2)*pow(delta, 3) + 8*pow(Delta, 2)*Delta_dd*pow(X, 2)*pow(delta, 4) + 10*pow(Delta, 2)*X*delta - 24*Delta*pow(Delta_d, 2)*pow(X, 2)*pow(delta, 4) + 8*Delta*Delta_d*X*pow(delta, 2) + 8*Delta*Delta_dd*X*pow(delta, 3) + 8*Delta - 18*pow(Delta_d, 2)*X*pow(delta, 3) + 12*Delta_d*delta + 2*Delta_dd*pow(delta, 2))/(16*pow(Delta, 4)*pow(X, 4)*pow(delta, 4) + 32*pow(Delta, 3)*pow(X, 3)*pow(delta, 3) + 24*pow(Delta, 2)*pow(X, 2)*pow(delta, 2) + 8*Delta*X*delta + 1);
    double dXdd_dDelta_d = 2*pow(X, 2)*(2*X*delta*(Delta + Delta_d*delta)*(Delta*X*delta + 1) + (2*Delta*X*delta + 1)*(2*Delta_d*X*pow(delta, 2) - 1))/pow(2*Delta*X*delta + 1, 3);
    double dXdd_dDelta_dd = -pow(X, 2)*delta/(2*Delta*X*delta + 1);

    return dXdd_ddelta + dXdd_dX*(dX_ddelta + dX_dDelta*Delta_d) + dXdd_dDelta*Delta_d + dXdd_dDelta_d*Delta_dd + dXdd_dDelta_dd*Delta_ddd;
}


double ResidualHelmholtzSAFTAssociating::g(double eta)
{
	return 0.5*(2-eta)/pow(1-eta,(int)3);
}    
double ResidualHelmholtzSAFTAssociating::dg_deta(double eta)
{
	return 0.5*(5-2*eta)/pow(1-eta,(int)4);
}
double ResidualHelmholtzSAFTAssociating::d2g_deta2(double eta)
{
    return 3*(3-eta)/pow(1-eta,(int)5);
}   
double ResidualHelmholtzSAFTAssociating::d3g_deta3(double eta)
{
	return 6*(7-2*eta)/pow(1-eta,(int)6);
}   
double ResidualHelmholtzSAFTAssociating::eta(double delta){
	return this->vbarn*delta;
}
double ResidualHelmholtzSAFTAssociating::A(double log_tau, double tau, double log_delta, double delta, int i)
{
	double X = this->X(delta, this->Deltabar(tau, delta));
    return this->m*this->a*((log(X)-X/2.0+0.5));
}
double ResidualHelmholtzSAFTAssociating::dA_dDelta(double log_tau, double tau, double log_delta, double delta, int i)
{
	double X = this->X(delta, this->Deltabar(tau, delta));
    return this->m*this->a*(1/X-0.5)*this->dX_ddelta(tau, delta);
}
double ResidualHelmholtzSAFTAssociating::dA_dTau(double log_tau, double tau, double log_delta, double delta, int i)
{
	double X = this->X(delta, this->Deltabar(tau, delta));
    return this->m*this->a*(1/X-0.5)*this->dX_dtau(tau, delta);
}
double ResidualHelmholtzSAFTAssociating::d2A_dTau2(double log_tau, double tau, double log_delta, double delta, int i)
{
	double X = this->X(delta, this->Deltabar(tau, delta));
	double X_tau = this->dX_dtau(tau, delta);
    double X_tautau = this->d2X_dtau2(tau, delta);
    return this->m*this->a*((1/X-0.5)*X_tautau-pow(X_tau/X, 2));
}
double ResidualHelmholtzSAFTAssociating::d2A_dDelta2(double log_tau, double tau, double log_delta, double delta, int i)
{
	double X = this->X(delta, this->Deltabar(tau, delta));
	double X_delta = this->dX_ddelta(tau, delta);
    double X_deltadelta = this->d2X_ddelta2(tau, delta);
    return this->m*this->a*((1/X-0.5)*X_deltadelta-pow(X_delta/X,2));
}
double ResidualHelmholtzSAFTAssociating::d2A_dDelta_dTau(double log_tau, double tau, double log_delta, double delta, int i)
{
	double X = this->X(delta, this->Deltabar(tau, delta));
	double X_delta = this->dX_ddelta(tau, delta);
    double X_deltadelta = this->d2X_ddelta2(tau, delta);
    double X_tau = this->dX_dtau(tau, delta);
    double X_deltatau = this->d2X_ddeltadtau(tau, delta);
    return this->m*this->a*((-X_tau/X/X)*X_delta+X_deltatau*(1/X-0.5));
}
double ResidualHelmholtzSAFTAssociating::d3A_dTau3(double log_tau, double tau, double log_delta, double delta, int i)
{
	double X = this->X(delta, this->Deltabar(tau, delta));
    double X_t = this->dX_dtau(tau, delta);
    double X_tt = this->d2X_dtau2(tau, delta);
    double X_ttt = this->d3X_dtau3(tau, delta);
    return this->m*this->a*((1/X-1.0/2.0)*X_ttt+(-X_t/pow(X,(int)2))*X_tt-2*(pow(X,(int)2)*(X_t*X_tt)-pow(X_t,(int)2)*(X*X_t))/pow(X,(int)4));
}
double ResidualHelmholtzSAFTAssociating::d3A_dDelta_dTau2(double log_tau, double tau, double log_delta, double delta, int i)
{
	double X = this->X(delta, this->Deltabar(tau, delta));
    double X_t = this->dX_dtau(tau, delta);
    double X_d = this->dX_ddelta(tau, delta);
    double X_tt = this->d2X_dtau2(tau, delta);
    double X_dt = this->d2X_ddeltadtau(tau, delta);
    double X_dtt = this->d3X_ddeltadtau2(tau, delta);
    return this->m*this->a*((1/X-1.0/2.0)*X_dtt-X_d/pow(X,(int)2)*X_tt-2*(pow(X,(int)2)*(X_t*X_dt)-pow(X_t,(int)2)*(X*X_d))/pow(X,(int)4));
}
double ResidualHelmholtzSAFTAssociating::d3A_dDelta2_dTau(double log_tau, double tau, double log_delta, double delta, int i)
{
	double X = this->X(delta, this->Deltabar(tau, delta));
    double X_t = this->dX_dtau(tau, delta);
    double X_d = this->dX_ddelta(tau, delta);
    double X_dd = this->d2X_ddelta2(tau, delta);
    double X_dt = this->d2X_ddeltadtau(tau, delta);
    double X_ddt = this->d3X_ddelta2dtau(tau, delta);
    return this->m*this->a*((1/X-1.0/2.0)*X_ddt-X_t/pow(X,(int)2)*X_dd-2*(pow(X,(int)2)*(X_d*X_dt)-pow(X_d,(int)2)*(X*X_t))/pow(X,(int)4));
}
double ResidualHelmholtzSAFTAssociating::d3A_dDelta3(double log_tau, double tau, double log_delta, double delta, int i)
{
	double X = this->X(delta, this->Deltabar(tau, delta));
    double X_d = this->dX_ddelta(tau, delta);
    double X_dd = this->d2X_ddelta2(tau, delta);
    double X_ddd = this->d3X_ddelta3(tau, delta);
    return this->m*this->a*((1/X-1.0/2.0)*X_ddd-X_d/pow(X,(int)2)*X_dd-2*(pow(X,(int)2)*(X_d*X_dd)-pow(X_d,(int)2)*(X*X_d))/pow(X,(int)4));
}