//
// Created by thomae on 20/06/2023.
//

#include "LiftingPolicySolver.h"

#include <utility>

namespace robust_model {

SingleDirectionBreakPoints::SingleDirectionBreakPoints(helpers::SmartIndex<SingleDirectionBreakPoints> const& id,
                                                       SingleDirectionBreakPoints::BreakPointsSeries break_points,
                                                       SingleDirectionBreakPoints::BreakPointDirection break_point_direction)
        : IndexedObject<SingleDirectionBreakPoints>(id),
          _break_points(std::move(break_points)),
          _break_point_direction(std::move(break_point_direction)) {}

bool SingleDirectionBreakPoints::simple_axis_aligned() const {
    if (break_point_direction().scaled_variables().size() != 1) {
        return false;
    }
    auto const& svar = break_point_direction().scaled_variables().front();
    return (svar.scale() == 1) and (svar.variable()->id().raw_id() == id().raw_id());
}

SingleDirectionBreakPoints::BreakPointsSeries const& SingleDirectionBreakPoints::break_points() const {
    return _break_points;
}

SingleDirectionBreakPoints::BreakPointDirection const& SingleDirectionBreakPoints::break_point_direction() const {
    return _break_point_direction;
}

UncertaintyVariable::Index const& SingleDirectionBreakPoints::axis_direction() const {
    helpers::exception_check(simple_axis_aligned(), "Can only get axis direction of aligned direction!");
    return break_point_direction().scaled_variables().front().variable()->id();
}

std::vector<double>
SingleDirectionBreakPoints::lifted_uncertainty_realization(
        UncertaintyRealization const& uncertainty_realization) const {
    double const directed_realization = break_point_direction().value(uncertainty_realization);
    std::vector<double> lifted_realization(break_points().size() + 1);
    for (size_t i = 0; i <= break_points().size(); ++i) {
        double const lb = previous_break_point(i);
        double const ub = break_point(i);
        lifted_realization.at(i) = std::max(std::min(directed_realization, ub) - lb, 0.);
    }
    return lifted_realization;
}

SingleDirectionBreakPoints::BreakPoint SingleDirectionBreakPoints::previous_break_point(size_t i) const {
    return (i > 0) ?
           break_points().at(i - 1) :
           break_point_direction().lb();
}

SingleDirectionBreakPoints::BreakPoint SingleDirectionBreakPoints::break_point(size_t i) const {
    return (i < break_points().size())
           ? break_points().at(i) :
           break_point_direction().ub();
}

LiftingPolicySolver::LiftingPolicySolver(ROModel const& original_model) :
    solvers::AROPolicySolverBase(original_model) {}

void LiftingPolicySolver::add_break_points(LiftingPolicySolver::BreakPointsSeries const& break_points,
                                           LiftingPolicySolver::BreakPointDirection const& break_point_direction) {
//    helpers::exception_check(not _breakpoint_tightening,
//                             "If domination breakpoints tightening shall be used, no additional breakpoints may be set!");
//TODO Add an exception check for validity of breakpoints - i.e. strictly increasing!
    auto const id = base_add_object(break_points, break_point_direction);
    _all_simple_axis_aligned = _all_simple_axis_aligned and id->simple_axis_aligned();
}

void LiftingPolicySolver::add_domination_motivated_break_points() {
    helpers::exception_check(num_objects() == 0,
                             "Only add domination break points, when no other breakpoints were set!");
    auto dm = DominatingPolicySolver(model(), robust_model::DominatingPolicySolver::Mode::THOMAE);
    dm.build_vertices();
    auto const base_vertex = dm.base_vertex();
    for (auto const& uvar: model().uncertainty_variables()) {
        if (base_vertex->vertex_value(uvar.id()) > uvar.lb() and base_vertex->vertex_value(uvar.id()) < uvar.ub()) {
            add_break_points(BreakPointsSeries{base_vertex->vertex_value(uvar.id())},
                             BreakPointDirection{uvar.reference()});
        } else {
            add_break_points(BreakPointsSeries{}, BreakPointDirection{uvar.reference()});
        }
    }
}

void LiftingPolicySolver::add_equidistant_breakpoints(size_t num_pieces) {
    for (auto const& uvar: model().uncertainty_variables()) {
        BreakPointsSeries breakpoints;
        for (size_t breakpoint_id = 1; breakpoint_id < num_pieces; ++breakpoint_id) {
            breakpoints.emplace_back(
                    uvar.lb() + double(breakpoint_id) * (uvar.ub() - uvar.lb()) / double(num_pieces));
        }
        add_break_points(breakpoints,
                         BreakPointDirection{uvar.reference()});
    }
}

void LiftingPolicySolver::build_implementation() {
    helpers::exception_check(not built(), "Only build model once!");
    if (_all_simple_axis_aligned) {
        build_axis_aligned_model();
    } else {
        helpers::exception_throw("Not implemented yet!");
    }

    _affine_model = std::make_unique<AffineAdjustablePolicySolver>(_lifted_model);
    affine_model().build();
}

void LiftingPolicySolver::build_axis_aligned_model() {
    helpers::exception_check(_all_simple_axis_aligned, "Not all break points are simple axis aligned!");

    axis_aligned_lifted_model_add_lifted_uncertainty_variables();
    axis_aligned_lifted_model_add_retracted_uncertainty_constraints();

    add_decision_variables();
    axis_aligned_lifted_model_add_retracted_constraints();
    axis_aligned_lifted_model_add_retracted_objective();

    if (all_rotational_invariant_axis_aligned_breakpoints() and
        model().uncertainty_set().rotational_invariant() and
        _breakpoint_tightening) {
        add_rotational_invariant_tightening_constraints();
    }
}


void LiftingPolicySolver::add_decision_variables() {
    for (auto const& dvar: model().decision_variables()) {
        _lifted_decision_variables.emplace_back(
                lifted_model().add_decision_variable(dvar.name(),
                                                    (dvar.has_period()) ? dvar.period() : std::optional<period_id>{},
                                                    dvar.lb(), dvar.ub()));
    }
}

void LiftingPolicySolver::axis_aligned_lifted_model_add_lifted_uncertainty_variables() {
    for (auto const& break_point_series: objects()) {
        _lifted_uncertainty_variables.emplace_back();
        auto const& udirection = break_point_series.axis_direction();
        for (size_t i = 0; i <= break_point_series.break_points().size(); ++i) {
            BreakPoint break_point = break_point_series.break_point(i);
            BreakPoint prev_break_point = break_point_series.previous_break_point(i);
            _lifted_uncertainty_variables.back().emplace_back(
                    lifted_model().add_uncertainty_variable(
                            udirection->name() + "_L" + std::to_string(i),
                            udirection->has_period() ? udirection->period() : std::optional<period_id>{},
                            0,
                            break_point - prev_break_point
                    )
            );
            if (i > 0) {
                lifted_model().add_uncertainty_constraint(
                        _lifted_uncertainty_variables.back().at(i) /
                        _lifted_uncertainty_variables.back().at(i).ub()
                        <=
                        _lifted_uncertainty_variables.back().at(i - 1) /
                        _lifted_uncertainty_variables.back().at(i - 1).ub(),
                        "BoundLiftedWithPrevious" + udirection->name() + "_" + std::to_string(i)
                );
            }
        }
        _lifted_uncertainty_retractions.emplace_back(
                udirection->lb(),
                LinearExpression<UncertaintyVariable::Reference>::sum(_lifted_uncertainty_variables.back())
        );
    }
}

void LiftingPolicySolver::axis_aligned_lifted_model_add_retracted_uncertainty_constraints() {
    for (auto const union_uncertainty_set: model().uncertainty_set().constraint_sets()) {
        if (union_uncertainty_set.raw_id() > 0) {
            lifted_model().add_uncertainty_constraint_set();
        }
        auto const lifted_union_uncertainty_set = lifted_model().uncertainty_set().constraint_sets().back();
        for (auto const& uconstr: model().uncertainty_set().uncertainty_constraints(union_uncertainty_set)) {
            lifted_model().add_uncertainty_constraint(
                    uconstr.substitute<UncertaintyVariable>(
                            _lifted_uncertainty_retractions
                    ),
                    lifted_union_uncertainty_set);
        }
    }
}

void LiftingPolicySolver::axis_aligned_lifted_model_add_retracted_constraints() {
    for (auto const& roconstr: model().constraints()) {
        lifted_model().add_constraint(
                {
                        roconstr.name(),
                        roconstr.sense(),
                        roconstr.expression().substitute(_lifted_decision_variables,
                                                         _lifted_uncertainty_retractions)
                });
    }
}

void LiftingPolicySolver::axis_aligned_lifted_model_add_retracted_objective() {
    auto const& old_obj = model().objective();
    lifted_model().set_objective(
            old_obj.expression().substitute(_lifted_decision_variables, _lifted_uncertainty_retractions),
            old_obj.sense());
}

void LiftingPolicySolver::add_rotational_invariant_tightening_constraints() {
    helpers::exception_check(all_rotational_invariant_axis_aligned_breakpoints() and
                             model().uncertainty_set().rotational_invariant(),
                             "Only add rotational invariant tightening for rotational invariant sets!");
    size_t const lift_vars_per_variable = _lifted_uncertainty_variables.front().size();

    if (all_symmetric_axis_aligned_breakpoints() and model().uncertainty_set().symmetric()) {
        AffineExpression<UncertaintyVariable::Reference> lhs;
        for (size_t i = 0; i < lift_vars_per_variable / 2; ++i) {
            double break_point = -objects().begin()->break_point(i);
            for (auto const& lifted_vars: _lifted_uncertainty_variables) {
                lhs += lifted_vars.at(i).ub() - lifted_vars.at(i);
                lhs += lifted_vars.at(lift_vars_per_variable - i - 1);
            }
            lifted_model().add_uncertainty_constraint(
                    lhs <= max_rotational_invariant_outside_budget(break_point),
                    "RotationalBound_BP" + std::to_string(i)
            );
        }
    }
    if (model().uncertainty_set().non_negative()) {
        LinearExpression<UncertaintyVariable::Reference> lhs;
        for (int i = int(lift_vars_per_variable) - 1; i >= 0; --i) {
            double break_point = objects().begin()->previous_break_point(i);
            for (auto const& lifted_vars: _lifted_uncertainty_variables) {
                lhs += lifted_vars.at(i);
            }
            lifted_model().add_uncertainty_constraint(
                    lhs <= max_rotational_invariant_outside_budget(break_point),
                    "RotationalBound_BP" + std::to_string(i)
            );
        }
    }
}

double LiftingPolicySolver::max_rotational_invariant_outside_budget(double bound) const {
    double max_outside_budget = 0;
    for (size_t k = 1; k <= model().num_uvars(); ++k) {
        max_outside_budget = std::max(max_outside_budget,
                                      model().uncertainty_set().max_one_norm_k_active(k) - double(k) * bound);
    }
    return max_outside_budget;
}

void LiftingPolicySolver::solve_implementation() {
    helpers::exception_check(built(), "Can only solve built model!");
    set_parameters_to_other(affine_model());
    affine_model().solve();
    set_results_from_other(affine_model());
}

SolutionRealization LiftingPolicySolver::specific_solution(std::vector<double> const& uncertainty_realization) const {
    return {model(),
            uncertainty_realization,
            affine_model().specific_solution(lifted_uncertainty_realization(uncertainty_realization)).solutions()};
}

void LiftingPolicySolver::add_average_lifted_vertex_reoptimization() {
    helpers::exception_check(built(), "Only add reoptimization objective, when base model is build!");
    auto dm = DominatingPolicySolver(model(), robust_model::DominatingPolicySolver::Mode::THOMAE);
    dm.build_vertices();
    double scale = dm.approximation_factor();
    std::vector<double> average_lifted_uncertainty(lifted_model().num_uvars(), 0.);
    for (size_t vertex = 0; vertex < dm.num_vertices(); ++vertex) {
        std::vector<double> scaled_vertex;
        for (double const v: dm.vertex_id(vertex)->vertex_vector()) {
            scaled_vertex.emplace_back(v / scale);
        }
        auto const lifted_vertex = lifted_uncertainty_realization(scaled_vertex);
        for (size_t i = 0; i < lifted_model().num_uvars(); ++i) {
            average_lifted_uncertainty.at(i) += lifted_vertex.at(i) / double(dm.num_vertices());
        }
    }
    affine_model().add_reoptimization_objective_for_realization(UncertaintyRealization(average_lifted_uncertainty));
}

std::vector<double>
LiftingPolicySolver::lifted_uncertainty_realization(std::vector<double> const& uncertainty_realization) const {
    UncertaintyRealization const realization(uncertainty_realization);
    std::vector<double> lifted_uncertainty(lifted_model().num_uvars());
    for (auto const& break_point_series: objects()) {
        auto const break_direction_lifting = break_point_series.lifted_uncertainty_realization(realization);
        for (size_t i = 0; i <= break_point_series.break_points().size(); ++i) {
            auto const uvar = _lifted_uncertainty_variables.at(break_point_series.id().raw_id()).at(i);
            lifted_uncertainty.at(uvar.raw_id()) = break_direction_lifting.at(i);
        }
    }
    return lifted_uncertainty;
}

bool LiftingPolicySolver::all_rotational_invariant_axis_aligned_breakpoints() {
    if (not _all_simple_axis_aligned) {
        return false;
    }
    if (num_objects() == 0) {
        return true;
    }
    auto const& series = objects().front().break_points();
    for (auto const& break_points: objects()) {
        if (break_points.break_points().size() != series.size()) {
            return false;
        }
        for (size_t i = 0; i < series.size(); ++i) {
            if (series.at(i) != break_points.break_point(i)) {
                return false;
            }
        }
    }
    return true;
}

bool LiftingPolicySolver::all_symmetric_axis_aligned_breakpoints() {
    if (not _all_simple_axis_aligned) {
        return false;
    }
    if (num_objects() == 0) {
        return true;
    }
    for (auto const& break_points: objects()) {
        for (size_t i = 0; i < break_points.break_points().size(); ++i) {
            if (break_points.break_point(break_points.break_points().size() - i - 1) != -break_points.break_point(i)) {
                return false;
            }
        }
    }
    return true;
}

ROModel& LiftingPolicySolver::lifted_model() {
    return _lifted_model;
}

ROModel const& LiftingPolicySolver::lifted_model() const {
    return _lifted_model;
}

AffineAdjustablePolicySolver& LiftingPolicySolver::affine_model() {
    return *_affine_model;
}

AffineAdjustablePolicySolver const& LiftingPolicySolver::affine_model() const {
    return *_affine_model;
}

}