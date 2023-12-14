//
// Created by simon on 16/06/2021.
//

#ifndef ROBUSTOPTIMIZATION_SOLUTIONREALIZATION_H
#define ROBUSTOPTIMIZATION_SOLUTIONREALIZATION_H

#include "../ROModel.h"

namespace robust_model {

class SolutionRealization {
public:
    SolutionRealization(ROModel const& ro_model,
                     std::vector<double>  realization,
                     std::vector<double>  solutions);

    double value(UncertaintyVariable::Index uvar) const;

    double value(DecisionVariable::Index dvar) const;

    std::vector<double> const& realization() const;

    std::vector<double> const& solutions() const;

    std::string realization_string() const;

    std::string solution_string() const;

    bool feasible(double tolerance = 1e-3) const;

    double objective_value() const;

private:
    ROModel const& _ro_model;
    std::vector<double> const _realization;
    std::vector<double> const _solutions;
};

}

#endif //ROBUSTOPTIMIZATION_SOLUTIONREALIZATION_H
