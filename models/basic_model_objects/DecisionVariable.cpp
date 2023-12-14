//
// Created by simon on 14/06/2021.
//

#include "DecisionVariable.h"

namespace robust_model{

DominatingSolution::DominatingSolution(std::vector<double>  solutions) : _values_per_vertex(std::move(solutions)){}

std::vector<double> const& DominatingSolution::values_per_vertex() const {
    return _values_per_vertex;
}

std::string DominatingSolution::to_string() const {
    std::string s;
    for(auto const val : values_per_vertex()){
        s += std::to_string(val) + " ";
    }
    return s;
}

AffineSolution::AffineSolution(double constant,
                               std::map<UncertaintyVariable::Index, double>  dependent_scales):
        _constant(constant),
        _dependent_scales(std::move(dependent_scales)) {}

double AffineSolution::constant() const {
    return _constant;
}

std::map<UncertaintyVariable::Index, double> const& AffineSolution::dependent_scales() const {
    return _dependent_scales;
}

std::string AffineSolution::to_string() const {
    std::string s(std::to_string(constant()));
    for(auto const& [var, value] : dependent_scales()){
        s += " + " + std::to_string(value) + var->name();
    }
    return s;
}

DecisionVariable::DecisionVariable(Index id, std::string const& name,
                                   std::optional<period_id> p, double lb,
                                   double ub)
        : VariableBase<DecisionVariable>(id, name, lb, ub), _period(p) {
}

void DecisionVariable::add_dependency(UncertaintyVariable::Index dependency) {
    _dependencies.emplace_back(dependency);
}

void DecisionVariable::add_dependencies(std::vector<UncertaintyVariable::Reference> const& dependencies) {
    _dependencies.insert(_dependencies.end(), dependencies.begin(), dependencies.end());
}

std::vector<UncertaintyVariable::Reference> const& DecisionVariable::dependencies() const {
    return _dependencies;
}

bool DecisionVariable::has_period() const {
    return _period.has_value();
}

period_id DecisionVariable::period() const {
    return _period.value();
}

AffineSolution const& DecisionVariable::affine_solution() const {
    return _affine_solution.value();
}

void DecisionVariable::set_affine_solution(AffineSolution const& affine_solution) {
    _affine_solution = affine_solution;
}

DominatingSolution const& DecisionVariable::dominating_solution() const {
    return _dominating_solution.value();
}

void DecisionVariable::set_dominating_solution(DominatingSolution const& dominating_solution) {
    _dominating_solution = dominating_solution;
}

}
