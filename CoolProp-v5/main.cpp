
//#include <vld.h> 

#include "Backends/REFPROPMixtureBackend.h"
#include "Backends/REFPROPBackend.h"
#include <time.h>
#include "AbstractState.h"
#include "DataStructures.h"
#include <cstdio>
using namespace CoolProp;

#include "rapidjson/rapidjson_include.h"
#include "Fluids\FluidLibrary.h"

int main()
{

    if (0)
    {
        struct element
        {
            double d,t,ld;
            int l;
        };
        double n[] = {0.0125335479355233,                        7.8957634722828,                        -8.7803203303561,                        0.31802509345418,                        -0.26145533859358,                        -0.0078199751687981,                        0.0088089493102134,                        -0.66856572307965,                        0.20433810950965,                        -6.621260503968699e-005,                        -0.19232721156002,                        -0.25709043003438,                        0.16074868486251,                        -0.040092828925807,                        3.9343422603254e-007,                        -7.5941377088144e-006,                        0.00056250979351888,                        -1.5608652257135e-005,                        1.1537996422951e-009,                        3.6582165144204e-007,                        -1.3251180074668e-012,                        -6.2639586912454e-010,                        -0.10793600908932,                        0.017611491008752,                        0.22132295167546,                        -0.40247669763528,                        0.58083399985759,                        0.0049969146990806,                        -0.031358700712549,                        -0.74315929710341,                        0.4780732991548,                        0.020527940895948,                        -0.13636435110343,                        0.014180634400617,                        0.008332650488071301,                        -0.029052336009585,                        0.038615085574206,                        -0.020393486513704,                        -0.0016554050063734,                        0.0019955571979541,                        0.00015870308324157,                        -1.638856834253e-005,                        0.043613615723811,                        0.034994005463765,                        -0.076788197844621,                        0.022446277332006,                        -6.2689710414685e-005,                        -5.5711118565645e-010,                        -0.19905718354408,                        0.31777497330738,                        -0.11841182425981  }; 
        double d[] = {1,                         1,                        1,                        2,                        2,                        3,                        4,                        1,                        1,                        1,                        2,                        2,                        3,                        4,                        4,                        5,                        7,                        9,                        10,                        11,                        13,                        15,                        1,                        2,                        2,                        2,                        3,                        4,                        4,                        4,                        5,                        6,                        6,                        7,                        9,                        9,                        9,                        9,                        9,                        10,                        10,                        12,                        3,                        4,                        4,                        5,                        14,                        3,                        6,                        6,                        6                    };                   
        double t[] = {-0.5,                        0.875,                        1,                        0.5,                        0.75,                        0.375,                        1,                        4,                        6,                        12,                        1,                        5,                        4,                        2,                        13,                        9,                        3,                        4,                        11,                        4,                        13,                        1,                        7,                        1,                        9,                        10,                        10,                        3,                        7,                        10,                        10,                        6,                        10,                        10,                        1,                        2,                        3,                        4,                        8,                        6,                        9,                        8,                        16,                        22,                        23,                        23,                        10,                        50,                        44,                        46,                        50                    };
        double l[] = {0,                        0,                        0,                        0,                        0,                        0,                        0,                        1,                        1,                        1,                        1,                        1,                        1,                        1,                        1,                        1,                        1,                        1,                        1,                        1,                        1,                        1,                        2,                        2,                        2,                        2,                        2,                        2,                        2,                        2,                        2,                        2,                        2,                        2,                        2,                        2,                        2,                        2,                        2,                        2,                        2,                        2,                        3,                        3,                        3,                        3,                        4,                        6,                        6,                        6,                        6                    };
        double summer = 0;
        std::vector<element> elements;
        for (std::size_t i = 0; i < 51; ++i)
        {
            element el;
            el.d = d[i];
            el.t = t[i];
            el.l = l[i];
            el.ld = (double)l[i];
            elements.push_back(el);
        }

        long N = 1000000;
        double t1 = clock();
        for (std::size_t ii = 0; ii < N; ++ii)
        {
            double delta = 1.3, tau = 0.7;
            double log_tau  = log(tau), log_delta = log(delta);
            for (int jj = 0; jj < 51; ++jj)
            {
                double di = d[jj], ti = t[jj], lid = l[jj];
                int li = (int)lid;
                double pow_delta_li;
                if (li > 0){
                    pow_delta_li = pow(delta, li);
                    summer += (di-lid*pow_delta_li)*exp(ti*log_tau+(di-1)*log_delta-pow_delta_li);
                }
                else{
                    summer += di*exp(ti*log_tau+(di-1)*log_delta);    
                }
            }
        }
        double t2 = clock();
        double elap = (t2-t1)/CLOCKS_PER_SEC/((double)N)*1e6;
        printf("%g %g\n",elap, summer);
        
    }
    if (0)
    {
        double rhomass = 1, T = 300;

        AbstractState *Mix = AbstractState::factory("CORE-n-Propane");
        Mix->update(DmassT_INPUTS, rhomass, T);
        double p1 = Mix->p();

        AbstractState *MixRP = AbstractState::factory("REFPROP-propane");
        MixRP->update(DmassT_INPUTS, rhomass, T);
        double p2 = MixRP->p();

        double rr =0 ;
    }
    if (1)
    {


        int N = 2;
        std::vector<double> z(N, 1.0/N);
        double rhomass = 1, T = 300;

        AbstractState *Mix = AbstractState::factory("CORE-Ethane|n-Propane");
        Mix->set_mole_fractions(z);
        Mix->update(DmassT_INPUTS, rhomass, T);
        double p1 = Mix->p();

        AbstractState *MixRP = AbstractState::factory("REFPROP-Ethane|Propane");
        MixRP->set_mole_fractions(z);
        MixRP->update(DmassT_INPUTS, rhomass, T);
        double p2 = MixRP->p();

        double rr = 0;
    }
    if(1)
    {
        time_t t1,t2;
        
        std::size_t N = 1000000;
        AbstractState *State = AbstractState::factory("CORE-Water");
        double p = State->p();
        double summer = 0;
        t1 = clock();
        for (std::size_t ii = 0; ii < N; ++ii)
        {
            //AbstractState *State = new REFPROPBackend("Methane");
            //summer += EOS->dalphar_dDelta(0.7,1.3);
            /*for (int i = 0; i < 50; i++)
            {
                summer += exp(1.3+1e-10*ii+1e-10*i);
            }
            summer += log(0.7-1e-10*ii);*/
            //summer += log(1.3);
            State->update(PT_INPUTS,101325,300);
            summer += State->p();
        }
        t2 = clock();
        delete State;
        double elap = ((double)(t2-t1))/CLOCKS_PER_SEC/((double)N)*1e6;
        printf("%g %g\n",elap, summer/((double)N));
        double eee = 0;
        return 0;

    }
    if (0)
    {
        AbstractState *State = AbstractState::factory("REFPROP-Methane|Ethane");
        std::vector<double> x(2,0.5);
        State->set_mole_fractions(x);
        State->update(DmassT_INPUTS,1,250);
        double hh = State->hmolar();
        double mu = State->viscosity();
        double sigma = State->surface_tension();
        delete State;
    }
    if (0)
    {
        time_t t1,t2;
        t1 = clock();
        long N = 100000;
        for (long ii = 0; ii < N; ii++)
        {
            AbstractState *State = AbstractState::factory("REFPROP-Methane");
            //AbstractState *State = new REFPROPBackend("Methane");
            delete State;
        }
        t2 = clock();
        double elap = ((double)(t2-t1))/CLOCKS_PER_SEC/((double)N)*1e6;
        printf("%g\n",elap);
    }

    if(0)
    {
        AbstractState *State = AbstractState::factory("REFPROP-Methane");
        State->update(DmassT_INPUTS,1,300);
        double hh = State->hmolar();
        double mu = State->viscosity();

        time_t t1,t2;
        t1 = clock();
        for (long ii = 0; ii < 1000000; ii++)
        {
            State->update(PQ_INPUTS,300000,1-ii*1e-6);
            //State->update(DmassT_INPUTS,1-ii*1e-10,180);
            //double hh1 = State->hmolar();
            //double mu2 = State->viscosity();
        }
        t2 = clock();
        double elap = ((double)(t2-t1))/CLOCKS_PER_SEC;
        printf("%g\n",elap);

        //double sigma = State->surface_tension();
        delete State;
    }	
}