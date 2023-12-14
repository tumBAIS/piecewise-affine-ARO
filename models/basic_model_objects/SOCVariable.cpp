//
// Created by simon on 14/06/2021.
//

#include "SOCVariable.h"

namespace robust_model {

DualInformation::DualInformation(DualInformation::DualType type, size_t constraint_number) :
        _type(type),
        _constraint_number(constraint_number) {

}

DualInformation::DualInformation(DualInformation::DualType type, size_t constraint_number, size_t normed_vector_entry) :
        _type(type),
        _constraint_number(constraint_number) {

}

DualInformation::DualInformation(DualInformation::DualType type, helpers::SmartIndex<SOCVariable> variable) :
        _type(type),
        _variable(variable) {

}

bool DualInformation::is_constraint_dual() const{
    return _type == DualType::CONSTRAINT;
}

bool DualInformation::is_normed_vector_dual() const{
    return _type == DualType::CONSTRAINT_NORMED_VECTOR;
}

bool DualInformation::is_lower_bound_dual() const{
    return _type == DualType::VARIABLE_LB;
}

bool DualInformation::is_upper_bound_dual() const{
    return _type == DualType::VARIABLE_UB;
}

size_t DualInformation::constraint_number() const{
    return _constraint_number.value();
}

size_t DualInformation::normed_vector_entry() const{
    return _normed_vector_entry.value();
}

helpers::SmartIndex<SOCVariable> DualInformation::bounded_dual_variable() const{
    return _variable.value();
}

SOCVariable::SOCVariable(Index id, std::string const& name, double lb, double ub, VariableType type) :
        SOCVariable(id, name, lb, ub, {}, type) {}

SOCVariable::SOCVariable(Index id, std::string const& name, double lb, double ub,
                         DualInformation const& dual_info) :
        SOCVariable(id, name, lb, ub, std::optional<DualInformation>{dual_info}, VariableType::Continuous) {}

SOCVariable::SOCVariable(Index id, std::string const& name, double lb, double ub,
                         std::optional<DualInformation> const& dual_info, VariableType type) :
        VariableBase(id, name, lb, ub), _type(type), _dual_info(dual_info) {}

void SOCVariable::set_solution(double sol) {
    _solution = sol;
}

void SOCVariable::invalidate_solution() {
    _solution = {};
}

bool SOCVariable::has_solution() const {
    return _solution.has_value();
}

double SOCVariable::solution() const {
    helpers::exception_check(_solution.has_value(), "No solution found yet!");
    return _solution.value();
}

VariableType const SOCVariable::type() const {
    return _type;
}

bool SOCVariable::has_dual_info() const {
    return _dual_info.has_value();
}

DualInformation const& SOCVariable::dual_information() const {
    return _dual_info.value();
}

}
