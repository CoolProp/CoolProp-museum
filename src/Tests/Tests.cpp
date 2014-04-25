
#include "Tests.h"

#include "AbstractState.h"
#include "DataStructures.h"

#include <assert.h>

#if defined ENABLE_CATCH

    #define CATCH_CONFIG_RUNNER
    #include "Catch.hpp"

    class ConsistencyFixture
    {
    protected:
        CoolProp::AbstractState *pState;
        int pair;
    public:
        ConsistencyFixture()
        {
            pState = NULL;
        }
        ~ConsistencyFixture()
        {
            delete pState;
        }
        void set_backend(std::string backend_string)
        {
            pState = CoolProp::AbstractState::factory(backend_string);
        }
        void set_pair(int pair)
        {
            this->pair = pair;
        }
        bool single_phase_consistency_check(long double T, long double p)
        {
            CoolProp::AbstractState &State = *pState;

            // Start with T,P as inputs, cycle through all the other pairs that are supported
            State.update(CoolProp::PT_INPUTS, p, T);
            long double hmolar = State.hmolar(), smolar = State.smolar(), rhomolar = State.rhomolar(), umolar = State.umolar();

            switch (pair)
            {
            /// In this group, T is one of the known inputs, iterate for the other one (easy)
            case CoolProp::HmolarT_INPUTS:
                State.update(pair, hmolar, T); break;
            case CoolProp::SmolarT_INPUTS:
                State.update(pair, smolar, T); break;
            case CoolProp::TUmolar_INPUTS:
                State.update(pair, umolar, T); break;
            case CoolProp::DmolarT_INPUTS:
                State.update(pair, rhomolar, T); break;

            /// In this group, p is one of the known inputs (a little less easy)
            case CoolProp::DmolarP_INPUTS:
                State.update(pair, rhomolar, p); break;
            case CoolProp::HmolarP_INPUTS:
                State.update(pair, hmolar, p); break;
            case CoolProp::PSmolar_INPUTS:
                State.update(pair, p, smolar); break;
            case CoolProp::PUmolar_INPUTS:
                State.update(pair, p, umolar); break;

            /// In this group, density is one of the known inputs (a bit harder)
            case CoolProp::DmolarHmolar_INPUTS:
                State.update(pair, rhomolar, hmolar); break;
            case CoolProp::DmolarSmolar_INPUTS:
                State.update(pair, rhomolar, smolar); break;
            case CoolProp::DmolarUmolar_INPUTS:
                State.update(pair, rhomolar, umolar); break;

            case CoolProp::HmolarSmolar_INPUTS:
                State.update(pair, hmolar, smolar); break;
            case CoolProp::SmolarUmolar_INPUTS:
                State.update(pair, smolar, umolar); break;
            }

            // Make sure we end up back at the same temperature and pressure we started out with
            if(fabs(T-State.T()) > 1e-2) throw CoolProp::ValueError(format("Error on T [%g K] is greater than 1e-2",fabs(State.T()-T)));
            if(fabs(p-State.p())/p > 1e-2)  throw CoolProp::ValueError(format("Error on p [%g %%] is greater than 1e-2",fabs(p-State.p())/p ));
            return true;
        }
    };

    TEST_CASE_METHOD(ConsistencyFixture, "Test all input pairs for water using all valid backends", "[]")
    {
        set_backend("REFPROP-Water");
        int inputs[] = {
            CoolProp::DmolarT_INPUTS, ///< Molar density in mol/m^3, Temperature in K
            CoolProp::HmolarT_INPUTS, ///< Enthalpy in J/mol, Temperature in K
            CoolProp::SmolarT_INPUTS, ///< Entropy in J/mol/K, Temperature in K
            CoolProp::TUmolar_INPUTS, ///< Internal energy in J/mol, Temperature in K             
            CoolProp::DmolarP_INPUTS, ///< Molar density in mol/m^3, Pressure in Pa
            CoolProp::HmolarP_INPUTS, ///< Enthalpy in J/mol, Pressure in Pa
            CoolProp::PSmolar_INPUTS, ///< Pressure in Pa, Entropy in J/mol/K 
            CoolProp::PUmolar_INPUTS, ///< Pressure in Pa, Internal energy in J/mol
            CoolProp::HmolarSmolar_INPUTS, ///< Enthalpy in J/mol, Entropy in J/mol/K
            CoolProp::SmolarUmolar_INPUTS, ///< Entropy in J/mol/K, Internal energy in J/mol
            CoolProp::DmolarHmolar_INPUTS, ///< Molar density in mol/m^3, Enthalpy in J/mol
            CoolProp::DmolarSmolar_INPUTS, ///< Molar density in mol/m^3, Entropy in J/mol/K
            CoolProp::DmolarUmolar_INPUTS, ///< Molar density in mol/m^3, Internal energy in J/mol
        };

        int N = sizeof(inputs)/sizeof(int);
        for (double T = 300; T < 400; T += 10)
        {
            for (int i = 0; i < N; i++)
            {
                int pair = inputs[i];
                std::string pair_desc = CoolProp::get_input_pair_short_desc(pair);
                set_pair(pair);
                CAPTURE(pair_desc);
                CAPTURE(T);
                CHECK_NOTHROW(single_phase_consistency_check(T, 10000));
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
        session.useConfigData(config);
        session.run();
    }

#endif