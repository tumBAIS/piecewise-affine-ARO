//
// Created by simon on 29/07/2022.
//

#ifndef ROBUSTOPTIMIZATION_SOSCONSTRAINT_H
#define ROBUSTOPTIMIZATION_SOSCONSTRAINT_H

#include "SOCVariable.h"

namespace robust_model {


class SOSConstraint {
public:
    explicit SOSConstraint(std::vector<SOCVariable::Reference> const& exclusive_variables);

    std::vector<SOCVariable::Reference> const& exclusive_variables() const;

    std::string to_string() const;

private:
    std::vector<SOCVariable::Reference> const _exclusive_variables;
};

}

#endif //ROBUSTOPTIMIZATION_SOSCONSTRAINT_H
