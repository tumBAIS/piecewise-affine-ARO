//
// Created by thomae on 12/8/2023.
//

#ifndef ROBUSTOPTIMIZATION_SOCSOLVERBASE_H
#define ROBUSTOPTIMIZATION_SOCSOLVERBASE_H

#include "../SolverBase.h"
#include "../../models/SOCModel.h"

namespace solvers {

class SOCSolverBase : public SolverBase{
public:
    explicit SOCSolverBase(robust_model::SOCModel& soc_model);

protected:
    robust_model::SOCModel const& soc_model() const;
    robust_model::SOCModel & non_const_soc_model();

private:
    robust_model::SOCModel& _soc_model;
};

}

#endif //ROBUSTOPTIMIZATION_SOCSOLVERBASE_H
