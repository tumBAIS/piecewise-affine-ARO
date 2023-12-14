//
// Created by thomae on 20/06/2023.
//

#ifndef ROBUSTOPTIMIZATION_LIFTINGPOLICYSOLVER_H
#define ROBUSTOPTIMIZATION_LIFTINGPOLICYSOLVER_H

#include "AffineAdjustablePolicySolver.h"
#include "DominatingPolicySolver.h"

namespace robust_model {

class SingleDirectionBreakPoints : public helpers::IndexedObject<SingleDirectionBreakPoints> {
public:
    using BreakPoint = double;
    using BreakPointsSeries = std::vector<BreakPoint>;
    using BreakPointDirection = LinearExpression<UncertaintyVariable::Reference>;

    SingleDirectionBreakPoints(helpers::SmartIndex<SingleDirectionBreakPoints> const& id,
                               BreakPointsSeries break_points, BreakPointDirection break_point_direction);

    bool simple_axis_aligned() const;

    BreakPointsSeries const& break_points() const;

    BreakPointDirection const& break_point_direction() const;

    UncertaintyVariable::Index const& axis_direction() const;

    BreakPoint previous_break_point(size_t i) const;

    BreakPoint break_point(size_t i) const;

    std::vector<double> lifted_uncertainty_realization(UncertaintyRealization const& uncertainty_realization) const;

private:
    BreakPointsSeries const _break_points;
    BreakPointDirection const _break_point_direction;
};

class LiftingPolicySolver :
        public solvers::AROPolicySolverBase,
        public helpers::IndexedObjectOwner<SingleDirectionBreakPoints> {
public:
    using BreakPoint = SingleDirectionBreakPoints::BreakPoint;
    using BreakPointsSeries = SingleDirectionBreakPoints::BreakPointsSeries;
    using BreakPointDirection = SingleDirectionBreakPoints::BreakPointDirection;

public:

    explicit LiftingPolicySolver(ROModel const& original_model);

    void add_break_points(BreakPointsSeries const& break_points, BreakPointDirection const& break_point_direction);

    void add_domination_motivated_break_points();

    void add_equidistant_breakpoints(size_t num_pieces);

    SolutionRealization specific_solution(std::vector<double> const& uncertainty_realization) const final;

    void add_average_lifted_vertex_reoptimization();

private:

    void solve_implementation() final;

    void build_implementation() final;

    void build_axis_aligned_model();

    void add_decision_variables();

    void axis_aligned_lifted_model_add_lifted_uncertainty_variables();

    void axis_aligned_lifted_model_add_retracted_uncertainty_constraints();

    void axis_aligned_lifted_model_add_retracted_constraints();

    void axis_aligned_lifted_model_add_retracted_objective();

    void add_rotational_invariant_tightening_constraints();

    double max_rotational_invariant_outside_budget(double bound) const;

    std::vector<double> lifted_uncertainty_realization(std::vector<double> const& uncertainty_realization) const;

    bool all_rotational_invariant_axis_aligned_breakpoints();

    bool all_symmetric_axis_aligned_breakpoints();

    ROModel & lifted_model();
    ROModel const& lifted_model() const;
    AffineAdjustablePolicySolver & affine_model();
    AffineAdjustablePolicySolver const& affine_model() const;

private:
    bool _all_simple_axis_aligned = true;
    bool _breakpoint_tightening = true;
    std::unique_ptr<AffineAdjustablePolicySolver> _affine_model;

    ROModel _lifted_model;
    std::vector<DecisionVariable::Reference> _lifted_decision_variables;
    std::vector<std::vector<UncertaintyVariable::Reference>> _lifted_uncertainty_variables;
    std::vector<AffineExpression<UncertaintyVariable::Reference>> _lifted_uncertainty_retractions;
};

}

#endif //ROBUSTOPTIMIZATION_LIFTINGPOLICYSOLVER_H
