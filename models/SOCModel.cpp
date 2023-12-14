//
// Created by simon on 20/05/2021.
//


#include "SOCModel.h"

#include <utility>
#include <cmath>
#include <algorithm>
#include "gurobi_c++.h"
#include "../solvers/soc_solvers/GurobiSOCSolver.h"

namespace robust_model {

SOCModel::SOCModel(std::string name) : _name(std::move(name)) {}

std::vector<SOCVariable> const& SOCModel::variables() const {
    return objects();
}

std::vector<SOCConstraint<SOCVariable>> const&
SOCModel::soc_constraints() const {
    return _soc_constraints;
}

std::vector<SOSConstraint> const& SOCModel::sos_constraints() const {
    return _sos_constraints;
}

SOCModel::Objective const&
SOCModel::objective(size_t const id) const {
    return _objectives.at(id);
}

std::vector<SOCModel::Objective> const& SOCModel::objectives() const {
    return _objectives;
}

void SOCModel::clear_and_set_objective(SOCModel::Objective const& objective) {
    _objectives.clear();
    add_objective(objective);
}

void
SOCModel::add_objective(Objective const& objective) {
    _objectives.emplace_back(objective);
}

bool SOCModel::is_multi_objective() const {
    return objectives().size() > 1;
}

bool SOCModel::is_continuous() const {
    for (auto const& var: variables()) {
        if (var.type() != VariableType::Continuous) {
            return false;
        }
    }
    return true;
}

bool SOCModel::all_affine() const {
    return std::all_of(soc_constraints().begin(), soc_constraints().end(),
                       [](auto c) { return c.expression().is_affine(); })
           and sos_constraints().empty();
}

SOCModel::SOCVariableReference
SOCModel::add_variable(std::string const& name, double lb, double ub, VariableType type) {
    return base_add_object(name, lb, ub, type)->reference();
}

SOCModel::SOCVariableReference
SOCModel::add_variable(std::string const& name, double lb, double ub, DualInformation const& info) {
    return base_add_object(name, lb, ub, info)->reference();
}

std::vector<SOCVariable::Reference>
SOCModel::add_variables(std::size_t n, std::string const& name, double lb, double ub, VariableType type) {
    std::vector<SOCVariable::Reference> new_variables;
    for (std::size_t i = 0; i < n; ++i) {
        new_variables.emplace_back(add_variable(name + "_" + std::to_string(i), lb, ub, type));
    }
    return new_variables;
}

void
SOCModel::add_constraint(SOCConstraint<SOCVariable> const& constraint) {
    _soc_constraints.emplace_back(constraint);
}

void SOCModel::add_sos_constraint(std::vector<SOCVariable::Reference> const& exclusive_variables) {
    _sos_constraints.emplace_back(exclusive_variables);
}

std::string const&
SOCModel::name() const {
    return _name;
}

std::string
SOCModel::full_string() const {
    std::string s = name();
    size_t obj_num = 0;
    for (auto const& obj: objectives()) {
        ++obj_num;
        s += "\nObjective " + std::to_string(obj_num) + ": " + obj.to_string();
    }
    for (auto const& constr: soc_constraints()) {
        s += "\n" + constr.to_string() +
             (constr.has_dual_value() ? (" (" + std::to_string(constr.dual_value()) + ")") : "");
    }
    for (auto const& constr: sos_constraints()) {
        s += "\n" + constr.to_string();
    }
    for (auto const& var: objects()) {
        s += "\n" + std::to_string(var.lb()) + " <= " + var.name() + " <= " + std::to_string(var.ub())
             + (var.has_solution() ? (" (" + std::to_string(var.solution()) + ")") : "");
    }
    return s;
}

void SOCModel::compute_dual() {
    invalidate_dual();
    _dual = std::make_unique<SOCModel>(name() + "Dual");
    AffineExpression<SOCVariable::Reference> dual_objective;
    std::vector<AffineExpression<SOCVariable::Reference>> dual_constraint_expressions(variables().size());
    for (size_t i = 0; i < soc_constraints().size(); ++i) {
        add_dual_of_expression(dual_objective, dual_constraint_expressions, i);
    }
    for (auto const& variable: variables()) {
        add_dual_of_var_bound(dual_objective, dual_constraint_expressions.at(variable.id().raw_id()), variable.id());
    }
    add_dual_of_objective(dual_objective, dual_constraint_expressions);
    for (auto const& var: variables()) {
        dual().add_constraint(dual_constraint_expressions[var.id().raw_id()] == 0,
                              "DVar_" + var.name());
    }
    if (objective().sense() == ObjectiveSense::MIN) {
        dual().add_objective({ObjectiveSense::MAX, dual_objective});
    } else {
        dual().add_objective({ObjectiveSense::MIN, -dual_objective});
    }
}

SOCModel const& SOCModel::dual() const {
    return *_dual;
}

SOCModel& SOCModel::dual() {
    return *_dual;
}

void SOCModel::invalidate_dual() {
    _dual.reset();
}

void SOCModel::add_dual_of_var_bound(AffineExpression<SOCVariable::Reference>& dual_objective,
                                     AffineExpression<SOCVariable::Reference>& dual_constraint_expression,
                                     SOCVariable::Index const& variable) {
    if (variable->lb() != NO_VARIABLE_LB) {
        auto const dv = dual().add_variable("DLB_" + variable->name(), NO_VARIABLE_LB, 0,
                                            DualInformation(DualInformation::DualType::VARIABLE_LB, variable));
        dual_objective += -variable->lb() * dv;
        dual_constraint_expression += 1 * dv;
    }
    if (variable->ub() != NO_VARIABLE_UB) {
        auto const dv = dual().add_variable("DUB_" + variable->name(), 0, NO_VARIABLE_UB,
                                            DualInformation(DualInformation::DualType::VARIABLE_UB, variable));
        dual_objective += -variable->ub() * dv;
        dual_constraint_expression += 1 * dv;
    }
}

void SOCModel::add_dual_of_expression(AffineExpression<SOCVariable::Reference>& dual_objective,
                                      std::vector<AffineExpression<SOCVariable::Reference>>& dual_constraint_expressions,
                                      size_t constraint_number) {
    auto const& constraint = soc_constraints().at(constraint_number);
    auto const dv = dual().add_variable("DE_" + constraint.name(), constraint.dual_lb(), constraint.dual_ub(),
                                        DualInformation(DualInformation::DualType::CONSTRAINT, constraint_number));
    dual_objective += constraint.soc_expression().affine().constant() * dv;
    for (auto const& svar: constraint.soc_expression().affine().linear().scaled_variables()) {
        dual_constraint_expressions[svar.variable().raw_id()] += svar.scale() * dv;
    }
    if (not constraint.soc_expression().is_affine()) {
        add_dual_of_normed_vector(dual_objective,
                                  dual_constraint_expressions, dv, constraint_number);
    }
}

void SOCModel::add_dual_of_normed_vector(AffineExpression<SOCVariable::Reference>& dual_objective,
                                         std::vector<AffineExpression<SOCVariable::Reference>>& dual_constraint_expressions,
                                         SOCVariable::Reference dv, size_t constraint_number) {
    auto const& constraint = soc_constraints().at(constraint_number);
    auto const& normed_vector = constraint.soc_expression().normed_vector();
    std::vector<SOCVariable::Reference> dus;
    for (size_t i = 0; i < normed_vector.normed_vector().size(); ++i) {
        auto const& affine = normed_vector.normed_vector().at(i);
        auto const du = dus.emplace_back(
                dual().add_variable("DNV_" + constraint.name() + "_" + std::to_string(i),
                                    NO_VARIABLE_LB, NO_VARIABLE_UB,
                                    DualInformation(DualInformation::DualType::CONSTRAINT_NORMED_VECTOR,
                                                    constraint_number, i)));
        dual_objective += affine.constant() * du;
        for (auto const& svar: affine.linear().scaled_variables()) {
            dual_constraint_expressions[svar.variable().raw_id()] += svar.scale() * du;
        }
    }
    dual().add_constraint(SOCExpression<SOCVariable>::norm(dus, dual_norm_type(normed_vector.norm_type())) <= dv,
                          "DSOC_" + constraint.name());
}

void SOCModel::add_dual_of_objective(AffineExpression<SOCVariable::Reference>& dual_objective,
                                     std::vector<AffineExpression<SOCVariable::Reference>>& dual_constraint_expressions) {
    double sign = (objective().sense() == ObjectiveSense::MIN) ? 1. : -1.;
    dual_objective += sign * objective().expression().constant();
    for (auto const& svar: objective().expression().linear().scaled_variables()) {
        dual_constraint_expressions.at(svar.variable().raw_id()) += sign * svar.scale();
    }
}

void SOCModel::set_dual_soc_constraint_values(std::vector<double> const& values) {
    helpers::exception_check(values.size() == _soc_constraints.size(),
                             "Need to set exactly one dual value for each soc constraint");
    for (size_t i=0; i<values.size(); ++i) {
        _soc_constraints.at(i).set_dual_value(values.at(i));
    }
}

void SOCModel::invalidate_solution() {
    for(auto & var : non_const_objects()){
        var.invalidate_solution();
    }
    for(auto & obj : _objectives){
        obj.invalidate_solution();
    }
    for(auto & constr : _soc_constraints){
        constr.invalidate_dual_value();
    }
}

}
