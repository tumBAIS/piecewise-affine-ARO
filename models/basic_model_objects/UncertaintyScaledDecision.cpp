//
// Created by simon on 09/05/2022.
//

#include "UncertaintyScaledDecision.h"


namespace robust_model {


UncertaintyScaledDecision::UncertaintyScaledDecision(UncertaintyVariable::Reference const& uncertainty_variable,
                                                     DecisionVariable::Reference const& decision_variable) :
        _uncertainty_variable(uncertainty_variable),
        _decision_variable(decision_variable) {}

UncertaintyScaledDecision::UncertaintyScaledDecision(DecisionVariable::Reference const& decision_variable,
                                                     UncertaintyVariable::Reference const& uncertainty_variable) :
        UncertaintyScaledDecision(uncertainty_variable, decision_variable) {}

UncertaintyVariable::Reference const& UncertaintyScaledDecision::uncertainty_variable() const {
    return _uncertainty_variable;
}

DecisionVariable::Reference const& UncertaintyScaledDecision::decision_variable() const {
    return _decision_variable;
}

std::string UncertaintyScaledDecision::to_string() const {
    return uncertainty_variable()->name() + "*" + decision_variable()->name();
}

double UncertaintyScaledDecision::lb() const {
    return std::min(
            {
                    _uncertainty_variable.lb() * _decision_variable.lb(),
                    _uncertainty_variable.lb() * _decision_variable.ub(),
                    _uncertainty_variable.ub() * _decision_variable.lb(),
                    _uncertainty_variable.ub() * _decision_variable.ub()
            }
    );
}

double UncertaintyScaledDecision::ub() const {
    return std::max(
            {
                    _uncertainty_variable.lb() * _decision_variable.lb(),
                    _uncertainty_variable.lb() * _decision_variable.ub(),
                    _uncertainty_variable.ub() * _decision_variable.lb(),
                    _uncertainty_variable.ub() * _decision_variable.ub()
            }
    );
}

}
