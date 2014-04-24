/*
 * AbstractBackend.h
 *
 *  Created on: 20 Dec 2013
 *      Author: jowr
 */

#ifndef HELMHOLTZEOSBACKEND_H_
#define HELMHOLTZEOSBACKEND_H_

#include <vector>
#include "HelmholtzEOSMixtureBackend.h"

namespace CoolProp {

class HelmholtzEOSBackend : public HelmholtzEOSMixtureBackend  {
private:
    /// DISABLED: Setting mole or mass fraction not allowed since pure and pseudo-pure have single component with mole fraction of 1
    void set_mole_fractions(const std::vector<double> &mole_fractions){throw NotImplementedError();};
    /// DISABLED: Setting mole or mass fraction not allowed since pure and pseudo-pure have single component with mole fraction of 1
    void set_mass_fractions(const std::vector<double> &mass_fractions){throw NotImplementedError();};
public:
    HelmholtzEOSBackend();
    HelmholtzEOSBackend(CoolPropFluid *pFluid){set_components(std::vector<CoolPropFluid*>(1,pFluid));};
    virtual ~HelmholtzEOSBackend(){};
};

} /* namespace CoolProp */
#endif /* HELMHOLTZEOSBACKEND_H_ */
