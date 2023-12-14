//
// Created by simon on 11.05.21.
//

#include "ROModel.h"
#include "../helpers/helpers.h"
#include "basic_model_objects/SolutionRealization.h"

#include <utility>
#include <algorithm>

namespace robust_model {

ROModel::ROModel(std::string name) : _name(std::move(name)), _uncertainty_set(*this) {}

std::vector<ROModel::DecisionReference>
ROModel::add_decision_variables(std::size_t n, std::string const& name, std::optional<period_id> p, double lb, double ub) {
    std::vector<DecisionReference> vars;
    for (std::size_t i = 0; i < n; ++i) {
        vars.emplace_back(add_decision_variable(name + std::to_string(i), p, lb, ub));
    }
    return vars;
}

std::vector<ROModel::DecisionReference>
ROModel::add_decision_variables_for_each_period(size_t n, std::string const& name, period_id first_period, double lb,
                                                double ub) {
    std::vector<DecisionReference> vars;
    for (std::size_t i = 0; i < n; ++i) {
        vars.emplace_back(add_decision_variable(name + std::to_string(i), first_period + period_id(i), lb, ub));
    }
    return vars;
}

ROModel::DecisionReference
ROModel::add_decision_variable(std::string const& name, std::optional<period_id> p, double lb, double ub) {
    auto const dvar = DecisionReference(base_add_object(name, p, lb, ub));
    if (not dvar->has_period()) {
        return dvar;
    }
    auto& non_const_dvar = object(dvar);
    for (auto const& uvar : uncertainty_variables()) {
        if (uvar.has_period() and uvar.period() <= dvar->period()) {
            non_const_dvar.add_dependency(uvar.id());
        }
    }
    return dvar;
}

ROModel::DecisionReference
ROModel::add_decision_variable(std::string const& name, std::vector<UncertaintyVariable::Reference> const& dependencies,
                               double lb, double ub) {
    auto const dvar = DecisionReference(base_add_object(name, std::optional<period_id>(), lb, ub));
    object(dvar).add_dependencies(dependencies);
    return dvar;
}

std::vector<ROModel::UncertaintyReference>
ROModel::add_uncertainty_variables(std::size_t n, std::string const& name, std::optional<period_id> p, double lb, double ub) {
    std::vector<UncertaintyReference> vars;
    for (std::size_t i = 0; i < n; ++i) {
        vars.emplace_back(add_uncertainty_variable(name + std::to_string(i), p, lb, ub));
    }
    return vars;
}

std::vector<ROModel::UncertaintyReference>
ROModel::add_uncertainty_variables_for_each_period(size_t n, std::string const& name, period_id first_period, double lb,
                                                   double ub) {
    std::vector<UncertaintyReference> vars;
    for (std::size_t i = 0; i < n; ++i) {
        vars.emplace_back(add_uncertainty_variable(name + std::to_string(i), first_period + period_id(i), lb, ub));
    }
    return vars;
}

ROModel::UncertaintyReference
ROModel::add_uncertainty_variable(std::string const& name, std::optional<period_id> p, double lb, double ub) {
    auto const uvar = _uncertainty_set.add_variable(name, p, lb, ub);
    if (not uvar->has_period()) {
        return uvar;
    }
    for (auto& dvar : non_const_objects()) {
        if (dvar.has_period() and uvar->period() <= dvar.period()) {
            dvar.add_dependency(uvar);
        }
    }
    return uvar;
}

void ROModel::add_constraint(RoConstraint const& constraint) {
    _constraints.emplace_back(constraint);
}

void ROModel::add_uncertainty_constraint(UncertaintySet::Constraint const& constraint) {
    _uncertainty_set.add_uncertainty_constraint(constraint);
}

UncertaintySetConstraintsSet::Index ROModel::add_uncertainty_constraint_set() {
    return _uncertainty_set.add_constraint_set();
}

void ROModel::add_uncertainty_constraint(UncertaintySet::Constraint const& constraint,
                                         UncertaintySetConstraintsSet::Index union_set) {
    _uncertainty_set.add_uncertainty_constraint(constraint, union_set);
}

void ROModel::set_objective(RoAffineExpression const& objective, ObjectiveSense sense) {
    _objective = {sense, objective};
}

std::string ROModel::full_string() const {
    std::string s = name();
    s += "\n" + objective().to_string();
    for (auto const& constr : constraints()) {
        s += "\n" + constr.to_string();
    }
    for (auto const& var:decision_variables()) {
        s += "\n" + std::to_string(var.lb()) + " <= " + var.name() + " p(" + std::to_string(var.period()) + ") <= " + std::to_string(var.ub());
    }
    s += "\n" + uncertainty_set().full_string();
    return s;
}

std::vector<DecisionVariable> const& ROModel::decision_variables() const {
    return objects();
}

std::vector<DecisionVariable::Index> const& ROModel::decision_variable_ids() const {
    return ids();
}

UncertaintySet const& ROModel::uncertainty_set() const {
    return _uncertainty_set;
}

void ROModel::add_special_uncertainty_constraint(UncertaintySet::SpecialSetType type, double budget) {
    _uncertainty_set.add_special_type_constraint(type, budget);
}

ROModel::Objective const& ROModel::objective() const {
    return _objective.value();
}

double ROModel::objective_value_for_solution(SolutionRealization const& solution) const {
    return objective().expression().value(solution);
}

std::vector<robust_model::ROModel::RoConstraint> const& ROModel::constraints() const {
    return _constraints;
}

std::string const& ROModel::name() const {
    return _name;
}

std::vector<UncertaintyVariable> const& ROModel::uncertainty_variables() const {
    return uncertainty_set().variables();
}

std::vector<UncertaintyVariable::Index> const& ROModel::uncertainty_variable_ids() const {
    return uncertainty_set().indices();
}

std::size_t ROModel::num_dvars() const {
    return num_objects();
}

std::size_t ROModel::num_uvars() const {
    return uncertainty_set().num_variables();
}

bool ROModel::solution_feasible(SolutionRealization const& solution, double tolerance) const {
    return std::all_of(constraints().begin(), constraints().end(),
                       [& solution, tolerance](RoConstraint const& c) {
                           return c.feasible(solution, tolerance);
                       });
}

}