
#include "Tests.h"

#include "AbstractState.h"
#include "DataStructures.h"

#include <assert.h>

#if defined ENABLE_CATCH

    #define CATCH_CONFIG_RUNNER
    #include "Catch.hpp"

    static int inputs[] = {
        CoolProp::DmolarT_INPUTS,
        /*CoolProp::SmolarT_INPUTS,
        CoolProp::HmolarT_INPUTS, 
        CoolProp::TUmolar_INPUTS, */
        
        CoolProp::DmolarP_INPUTS, 
        CoolProp::DmolarHmolar_INPUTS, 
        CoolProp::DmolarSmolar_INPUTS, 
        CoolProp::DmolarUmolar_INPUTS,
        
        /*
        CoolProp::HmolarP_INPUTS,
        CoolProp::PSmolar_INPUTS,
        CoolProp::PUmolar_INPUTS, 
        */

        /*
        CoolProp::HmolarSmolar_INPUTS, 
        CoolProp::HmolarUmolar_INPUTS, 
        CoolProp::SmolarUmolar_INPUTS 
        */
    };

    class ConsistencyFixture
    {
    protected:
        CoolProp::AbstractState *pState;
        int pair;
    public:
        ConsistencyFixture(){
            pState = NULL;
        }
        ~ConsistencyFixture(){
            delete pState;
        }
        void set_backend(std::string backend, std::string fluid_name){
            pState = CoolProp::AbstractState::factory(backend, fluid_name);
        }
        void set_pair(int pair){ 
            this->pair = pair;
        }
        bool single_phase_consistency_check(long double T, long double p)
        {
            CoolProp::AbstractState &State = *pState;

            // Start with T,P as inputs, cycle through all the other pairs that are supported
            State.update(CoolProp::PT_INPUTS, p, T);
            long double hmolar = State.hmolar(), smolar = State.smolar(), rhomolar = State.rhomolar(), umolar = State.umolar(), x1, x2;
            
            switch (pair)
            {
            /// In this group, T is one of the known inputs, iterate for the other one (easy)
            case CoolProp::HmolarT_INPUTS:
                x1 = hmolar; x2 = T;  break;
            case CoolProp::SmolarT_INPUTS:
                x1 = smolar; x2 = T; break;
            case CoolProp::TUmolar_INPUTS:
                x1 = T; x2 = umolar; break;
            case CoolProp::DmolarT_INPUTS:
                x1 = rhomolar; x2 = T; break;

            /// In this group, D is one of the known inputs, iterate for the other one (a little bit harder)
            case CoolProp::DmolarHmolar_INPUTS:
                x1 = rhomolar; x2 = hmolar; break;
            case CoolProp::DmolarSmolar_INPUTS:
                x1 = rhomolar; x2 = smolar; break;
            case CoolProp::DmolarUmolar_INPUTS:
                x1 = rhomolar; x2 = umolar; break;
            case CoolProp::DmolarP_INPUTS:
                x1 = rhomolar; x2 = p; break;

            /// In this group, p is one of the known inputs (a little less easy)
            case CoolProp::HmolarP_INPUTS:
                x1 = hmolar; x2 = p; break;
            case CoolProp::PSmolar_INPUTS:
                x1 = p; x2 = smolar; break;
            case CoolProp::PUmolar_INPUTS:
                x1 = p; x2 = umolar; break;

            case CoolProp::HmolarSmolar_INPUTS:
                x1 = hmolar; x2 = smolar; break;
            case CoolProp::SmolarUmolar_INPUTS:
                x1 = smolar; x2 = umolar; break;
            }
            CAPTURE(x1);
            CAPTURE(x2);
            //std::cout << format("input values were : %g, %g\n", x1, x2);
            State.update(pair, x1, x2);

            // Make sure we end up back at the same temperature and pressure we started out with
            if(fabs(T-State.T()) > 1e-2) throw CoolProp::ValueError(format("Error on T [%g K] is greater than 1e-2",fabs(State.T()-T)));
            if(fabs(p-State.p())/p*100 > 1e-2)  throw CoolProp::ValueError(format("Error on p [%g %%] is greater than 1e-2 %%",fabs(p-State.p())/p ));
            return true;
        }
    };

    TEST_CASE_METHOD(ConsistencyFixture, "Test all input pairs for CO2 using all valid backends", "[]")
    {
        set_backend("HEOS", "CO2");
        
        int N = sizeof(inputs)/sizeof(inputs[0]);
        for (double p = 600000; p < 800000000.0; p *= 5)
        {
            for (double T = 220; T < pState->Tmax(); T += 10)
            {
                for (int i = 0; i < N; i++)
                {
                    int pair = inputs[i];
                    std::string pair_desc = CoolProp::get_input_pair_short_desc(pair);
                    set_pair(pair);
                    CAPTURE(pair_desc);
                    CAPTURE(T);
                    CAPTURE(p);
                    CHECK_NOTHROW(single_phase_consistency_check(T, p));
                }
            }
        }
    }

    static Catch::Session session; // There must be exactly one instance

    int run_fast_tests()
    {
        Catch::ConfigData &config = session.configData();
        config.testsOrTags.clear();
        config.testsOrTags.push_back("[fast]");
        session.useConfigData(config);
        return session.run();
    }

    int run_not_slow_tests()
    {
        Catch::ConfigData &config = session.configData();
        config.testsOrTags.clear();
        config.testsOrTags.push_back("~[slow]");
        session.useConfigData(config);

        time_t t1, t2;
        t1 = clock();
        session.run();
        t2 = clock();
        printf("Elapsed time for not slow tests: %g s",(double)(t2-t1)/CLOCKS_PER_SEC);

        return 1;
    }

    int run_user_defined_tests(const std::vector<std::string> & tests_or_tags)
    {
        Catch::ConfigData &config = session.configData();
        config.testsOrTags.clear();
        for (unsigned int i = 0; i < tests_or_tags.size(); i++)
        {
            config.testsOrTags.push_back(tests_or_tags[i]);
        }
        session.useConfigData(config);

        time_t t1, t2;
        t1 = clock();
        session.run();
        t2 = clock();
        printf("Elapsed time for user defined tests: %g s",(double)(t2-t1)/CLOCKS_PER_SEC);

        return 1;
    }

    void run_tests()
    {
        Catch::ConfigData &config = session.configData();
        config.testsOrTags.clear();
        //config.shouldDebugBreak = true;
        session.useConfigData(config);
        session.run();
    }

#endif