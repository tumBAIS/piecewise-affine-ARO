//
// Created by thomae on 12/6/2023.
//

#ifndef ROBUSTOPTIMIZATION_AROPOLICYSOLVERBASE_H
#define ROBUSTOPTIMIZATION_AROPOLICYSOLVERBASE_H

#include "../../models/ROModel.h"
#include "../SolverBase.h"
#include "../../models/basic_model_objects/SolutionRealization.h"

namespace solvers {

class AROPolicySolverBase : public solvers::SolverBase {
public:
    explicit AROPolicySolverBase(robust_model::ROModel const& model);
    virtual robust_model::SolutionRealization specific_solution(std::vector<double> const& uncertainty_realization) const = 0;

protected:
    robust_model::ROModel const& model() const;

private:
    robust_model::ROModel const& _model;
};

}

#endif //ROBUSTOPTIMIZATION_AROPOLICYSOLVERBASE_H
