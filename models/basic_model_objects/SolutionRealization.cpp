//
// Created by simon on 16/06/2021.
//

#include "SolutionRealization.h"
#include "../ROModel.h"
#include <utility>

namespace robust_model{

SolutionRealization::SolutionRealization(ROModel const& ro_model, std::vector<double>  realization,
                                   std::vector<double>  solutions) : _ro_model(ro_model),
                                                                           _realization(std::move(realization)),
                                                                           _solutions(std::move(solutions)) {}

double SolutionRealization::value(UncertaintyVariable::Index uvar) const {
    return _realization.at(uvar.raw_id());
}

double SolutionRealization::value(DecisionVariable::Index dvar) const {
    return _solutions.at(dvar.raw_id());
}

std::vector<double> const& SolutionRealization::realization() const {
    return _realization;
}

std::vector<double> const& SolutionRealization::solutions() const {
    return _solutions;
}

std::string SolutionRealization::realization_string() const {
    std::string s;
    for(auto const& uvar : _ro_model.uncertainty_variables()){
        s += uvar.name() + ": " + std::to_string(value(uvar.id())) + "\n";
    }
    return s;
}

std::string SolutionRealization::solution_string() const {
    std::string s;
    for (auto const& dvar : _ro_model.decision_variables()) {
        s += dvar.name() + ": " + std::to_string(value(dvar.id())) + "\n";
    }
    return s;
}

bool SolutionRealization::feasible(double tolerance) const {
    return _ro_model.solution_feasible(*this, tolerance);
}

double SolutionRealization::objective_value() const {
    return _ro_model.objective_value_for_solution(*this);
}


}