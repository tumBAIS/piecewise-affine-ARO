//
// Created by simon on 29/07/2022.
//

#include "SOSConstraint.h"


namespace robust_model {


SOSConstraint::SOSConstraint(std::vector<SOCVariable::Reference> const& exclusive_variables) : _exclusive_variables(
        exclusive_variables) {}

std::vector<SOCVariable::Reference> const& SOSConstraint::exclusive_variables() const {
    return _exclusive_variables;
}

std::string SOSConstraint::to_string() const {
    std::string s = "SOS(";
    for(auto const& var:exclusive_variables()){
        s+=var->name() + ", ";
    }
    s.erase(s.size() - 2, 2);
    s += ")";
    return s;
}
}