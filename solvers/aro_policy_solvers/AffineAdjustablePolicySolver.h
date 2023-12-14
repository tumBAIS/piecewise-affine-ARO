//
// Created by simon on 21/05/2021.
//

#ifndef ROBUSTOPTIMIZATION_AFFINEADJUSTABLEPOLICYSOLVER_H
#define ROBUSTOPTIMIZATION_AFFINEADJUSTABLEPOLICYSOLVER_H

#include "../../models/ROModel.h"
#include "../../models/SOCModel.h"
#include "../../helpers/helpers.h"
#include "AROPolicySolverBase.h"
#include "../soc_solvers/GurobiSOCSolver.h"

namespace robust_model {

class AffineAdjustablePolicySolver : public solvers::AROPolicySolverBase {
public:
    explicit AffineAdjustablePolicySolver(ROModel const& model);

    AffineSolution affine_solution(DecisionVariable::Index dvar) const;

    SolutionRealization specific_solution(std::vector<double> const& uncertainty_realization) const final;

    void add_average_reoptimization();

    void add_reoptimization_objective_for_realization(UncertaintyRealization const& realization);

private:
    void solve_implementation() final;

    void build_implementation() final;

    void build_variables();

    void build_constraints();

    void build_objective();

    void add_ro_constraint(ROModel::RoConstraint const& ro_constr);

    AffineExpression<SOCVariable::Reference> add_rc_constraints_for_minimization(RoAffineExpression const& expr,
                                                                                 std::string const& name_addendum);

    void add_rc_constraints_for_equality(RoAffineExpression const& expr,
                                         std::string const& name_addendum);

    void
    add_dual_of_uncertainty_set(AffineExpression<SOCVariable::Reference>& dual_objective,
                                std::vector<AffineExpression<SOCVariable::Reference>>& dual_constraint_expressions,
                                UncertaintySetConstraintsSet::Index union_set_constraints,
                                std::string const& name_addendum);

    void
    add_dual_of_normed_vector(NormedAffineVector<UncertaintyVariable> const& normed_vector,
                              SOCVariable::Reference dv,
                              AffineExpression<SOCVariable::Reference>& dual_objective,
                              std::vector<AffineExpression<SOCVariable::Reference>>& dual_constraint_expressions,
                              std::string const& name_addendum);

    void
    add_dual_of_uvar_bounds(AffineExpression<SOCVariable::Reference>& dual_objective,
                            std::vector<AffineExpression<SOCVariable::Reference>>& dual_constraint_expressions,
                            std::string const& name_addendum);

    void
    add_dual_of_expression(AffineExpression<SOCVariable::Reference>& dual_objective,
                           std::vector<AffineExpression<SOCVariable::Reference>>& dual_constraint_expressions,
                           RoAffineExpression const& expr);


    SOCModel& soc_model();

    solvers::SOCSolverBase& soc_solver();

    std::vector<SOCVariable::Reference> const& adjustable_factors(DecisionVariable::Index id) const;

    SOCVariable::Reference const& adjustable_constant(DecisionVariable::Index id) const;


private:
    SOCModel _soc_model;
    std::unique_ptr<solvers::GurobiSOCSolver> _soc_solver;
    std::vector<std::vector<SOCVariable::Reference>> _adjustable_factors;
    std::vector<SOCVariable::Reference> _adjustable_constants;

};

}

#endif //ROBUSTOPTIMIZATION_AFFINEADJUSTABLEPOLICYSOLVER_H
