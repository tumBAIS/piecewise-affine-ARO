//
// Created by thomae on 12/8/2023.
//

#ifndef ROBUSTOPTIMIZATION_GUROBISOCSOLVER_H
#define ROBUSTOPTIMIZATION_GUROBISOCSOLVER_H

#include "gurobi_c++.h"
#include "SOCSolverBase.h"


namespace robust_model{
    class SOCModel;
}

namespace solvers {

//Constraint indices only correspond, when all constraints are affine!!!
//Similarly, the GRB model might have more variables, when not all constraints are affine!
class GurobiSOCSolver : public SOCSolverBase{
public:
    explicit GurobiSOCSolver(robust_model::SOCModel& soc_model);

    double value(robust_model::SOCVariable::Index const& id) const;

    void objectives_reset();

private:
    void solve_implementation() final;

    void build_implementation() final;

    void update_implementation() final;

    GRBModel & gurobi_model();

    void update_variables();
    void update_soc_constraints();
    void update_sos_constraints();
    void update_objectives();

    GRBLinExpr to_gurobi_linear(robust_model::AffineExpression<robust_model::SOCVariable::Reference> const& affine) const;

    helpers::VectorSlice<robust_model::SOCVariable> variables_to_add_grb() const;
    helpers::VectorSlice<robust_model::SOCConstraint<robust_model::SOCVariable>> soc_constraints_to_add_grb() const;
    helpers::VectorSlice<robust_model::SOSConstraint> sos_constraints_to_add_grb() const;
    helpers::VectorSlice<robust_model::ObjectiveBase<robust_model::AffineExpression<robust_model::SOCVariable::Reference>>> objectives_to_add_grb() const;

    void transfer_solution_to_soc_model();

private:

    std::unique_ptr<GRBEnv> _grb_env;
    std::unique_ptr<GRBModel> _grb_model;
    std::unique_ptr<std::vector<GRBVar>> _grb_vars;

    size_t _grb_next_var_to_add = 0;
    size_t _grb_next_constr_to_add = 0;
    size_t _grb_next_sos_constr_to_add = 0;
    size_t _grb_next_obj_to_add = 0;
};

}

#endif //ROBUSTOPTIMIZATION_GUROBISOCSOLVER_H
