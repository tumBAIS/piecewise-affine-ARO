//
// Created by thomae on 12/8/2023.
//

#include "SOCSolverBase.h"

namespace solvers{

robust_model::SOCModel const& SOCSolverBase::soc_model() const {
    return _soc_model;
}

robust_model::SOCModel& SOCSolverBase::non_const_soc_model() {
    return _soc_model;
}

SOCSolverBase::SOCSolverBase(robust_model::SOCModel& soc_model) : _soc_model(soc_model) {}

}
