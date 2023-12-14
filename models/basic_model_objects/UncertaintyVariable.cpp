//
// Created by simon on 14/06/2021.
//

#include "UncertaintyVariable.h"

namespace robust_model {

UncertaintyVariable::UncertaintyVariable(Index id,
                                         std::string const& name,
                                         std::optional<period_id> p, double lb, double ub) :
        VariableBase<UncertaintyVariable>(id, name, lb, ub), _period(p) {}

bool UncertaintyVariable::has_period() const {
    return _period.has_value();
}

period_id UncertaintyVariable::period() const {
    return _period.value();
}

bool UncertaintyVariable::is_0_1_var() const {
    return (lb()== 0.) and (ub() == 1.);
}

}
