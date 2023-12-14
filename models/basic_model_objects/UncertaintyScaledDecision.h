//
// Created by simon on 09/05/2022.
//

#ifndef ROBUSTOPTIMIZATION_UNCERTAINTYSCALEDDECISION_H
#define ROBUSTOPTIMIZATION_UNCERTAINTYSCALEDDECISION_H

#include "UncertaintyVariable.h"
#include "DecisionVariable.h"

namespace robust_model {

class SolutionRealization;

class UncertaintyScaledDecision {
public:
    UncertaintyScaledDecision(UncertaintyVariable::Reference const& uncertainty_variable,
                              DecisionVariable::Reference const& decision_variable);

    UncertaintyScaledDecision(DecisionVariable::Reference const& decision_variable,
                              UncertaintyVariable::Reference const& uncertainty_variable);

    UncertaintyVariable::Reference const& uncertainty_variable() const;

    DecisionVariable::Reference const& decision_variable() const;

    std::string to_string() const;

    template<class S>
    double value(S const& realization) const {
        return uncertainty_variable().value(realization) * decision_variable().value(realization);
    }

    double lb() const;
    double ub() const;

private:
    UncertaintyVariable::Reference _uncertainty_variable;
    DecisionVariable::Reference _decision_variable;
};

}

#endif //ROBUSTOPTIMIZATION_UNCERTAINTYSCALEDDECISION_H
