//
// Created by simon on 12.05.21.
//

#include "DominatingPolicySolver.h"

#include <utility>
#include "cmath"
#include "tuple"
#include "algorithm"
#include "../../helpers/helpers.h"

namespace robust_model {

Vertex::Vertex(Vertex::Index id, DominatingPolicySolver const& dominating_model, std::vector<double> vertex) :
        helpers::IndexedObject<Vertex>(id), _dominating_model(dominating_model), _vertex(std::move(vertex)) {

}

double Vertex::vertex_value(UncertaintyVariable::Index const& var) const {
    return _vertex.at(var.raw_id());
}

void Vertex::set_vertex_value(UncertaintyVariable::Index const& var, double const value) {
    _vertex.at(var.raw_id()) = value;
}

std::vector<double> const& Vertex::vertex_vector() const {
    return _vertex;
}

void Vertex::shift_vertex(UncertaintyVariable::Index const uvar, double const shift) {
    _vertex.at(uvar.raw_id()) += shift + -shift * vertex_value(uvar);
}

std::vector<UncertaintyVariable::Index> const& Vertex::dominated_uncertainties() const {
    return _dominated_uncertainties;
}

Vertex::Vertex(Index const& id, DominatingPolicySolver const& dominating_model, std::vector<double> const& vertex,
               std::vector<UncertaintyVariable::Index> const& dominated_uncertainties) :
        IndexedObject(id),
        _dominating_model(dominating_model),
        _vertex(vertex),
        _dominated_uncertainties(dominated_uncertainties) {}

std::string Vertex::to_string() const {
    std::string s;
    for (auto const& uvar: _dominating_model.model().uncertainty_variables()) {
        s += uvar.name() + ": " + std::to_string(vertex_value(uvar.id())) + "\n";
    }
    return s;
}

DominatingPolicySolver::DominatingPolicySolver(ROModel const& model, Mode mode) :
        solvers::AROPolicySolverBase(model),
        _mode(mode),
        _soc_model("Dominating " + model.name()) {
}

void DominatingPolicySolver::build_implementation() {
    helpers::exception_check(not built(), "Dominating model can only be built once!");
    auto const failed_assumption = first_failed_assumption();
    helpers::exception_check(failed_assumption.empty(), "Dominating Models only work for models of the form\n"
                                                        "MIN cTx \ns.t. Ax>=Bh+b,\n"
                                                        "where B, b are non negative and h >=0!\n"
                                                        "Failed Assumption: " + failed_assumption);
    if (num_vertices() == 0) {
        build_vertices();
    }
    build_variables();
    build_fix_period_constraints();
    build_constraints();
    build_objective();
    _soc_solver = std::make_unique<solvers::GurobiSOCSolver>(soc_model());
}

void DominatingPolicySolver::solve_implementation() {
    set_parameters_to_other(soc_solver());
    soc_solver().solve();
    set_results_from_other(soc_solver());
}

DominatingSolution DominatingPolicySolver::dominating_solution(DecisionVariable::Index const dvar) const {
    std::vector<double> solutions;
    for (auto const& v: VertexOwner::ids()) {
        solutions.emplace_back(dominating_variable(dvar, v)->solution());
    }
    return DominatingSolution(solutions);
}

SolutionRealization
DominatingPolicySolver::specific_solution(std::vector<double> const& uncertainty_realization) const {
    std::vector<double> solutions(model().num_dvars(), 0);
    double leftover_scale = 1;
    for (auto const v: VertexOwner::ids()) {
        if (v.raw_id() == base_vertex().raw_id())
            continue;
        double scale = vertex_scale(v, uncertainty_realization);
        leftover_scale -= scale;
        for (auto const& dvar: model().decision_variables()) {
            solutions.at(dvar.id().raw_id()) += scale * dominating_variable(dvar.id(), v)->solution();
        }
    }
//    helpers::exception_check(leftover_scale >= 0, "Total scales exceeded 1 (" +
//                                                  std::to_string(1 - leftover_scale) +
//                                                  ")! Did you use a point from within the uncertainty set?");
    for (auto const& dvar: model().decision_variables()) {
        solutions.at(dvar.id().raw_id()) +=
                leftover_scale * dominating_variable(dvar.id(), base_vertex())->solution();
    }
    return SolutionRealization(model(), uncertainty_realization, solutions);
}

double DominatingPolicySolver::vertex_scale(DominatingPolicySolver::VertexId vertex,
                                            std::vector<double> const& uncertainty_realization) const {
    auto const& dominated_uncertainties = DominatingPolicySolver::vertex(vertex).dominated_uncertainties();
    helpers::exception_check(not dominated_uncertainties.empty(),
                             "Can only find the scale of vertices that dominate at least one uncertainty.");
    auto const scale = [&](UncertaintyVariable::Index uvar) {
        double scale_nom = std::max(uncertainty_realization.at(uvar.raw_id())
                                    -
                                    vertex_solution_value(uvar, base_vertex()) *
                                    approach_scale(),
                                    0.);
        if (scale_nom == 0) {
            return 0.;
        }
        double scale_denom = (vertex_solution_value(uvar, vertex)
                              - vertex_solution_value(uvar, base_vertex()));
        return scale_nom / scale_denom;
    };
    return scale(
            *std::max_element(dominated_uncertainties.begin(), dominated_uncertainties.end(),
                              [scale](auto const& a, auto const& b) { return scale(a) < scale(b); })
    );

}

double DominatingPolicySolver::approach_scale() const {
    switch (mode()) {
        case Mode::BENTAL:
            switch (model().uncertainty_set().special_type()) {
                case UncertaintySet::SpecialSetType::BOX:
                    return 1;
                case UncertaintySet::SpecialSetType::BALL:
                    return .5;
                case UncertaintySet::SpecialSetType::BUDGET:
                    return 1;
                default:
                    helpers::exception_throw("Not implemented!");
            }
        default:
            return 1;
    }
}

SolutionRealization DominatingPolicySolver::specific_vertex_solution(std::size_t vertex) const {
    return specific_solution(DominatingPolicySolver::vertex(vertex_id(vertex)).vertex_vector());
}

void DominatingPolicySolver::shift_vertices(UncertaintyVariable::Index const uvar, double const shift) {
    helpers::exception_check(shift >= 0 and shift <= 1, "Can only shift between 0 and 1!");
    for (auto& v: VertexOwner::non_const_objects()) {
        v.shift_vertex(uvar, shift);
    }
}

DominatingPolicySolver::VertexId DominatingPolicySolver::vertex_id(std::size_t id) const {
    return VertexId(id, this);
}

Vertex const& DominatingPolicySolver::vertex(DominatingPolicySolver::VertexId id) const {
    return VertexOwner::object(id);
}

Vertex& DominatingPolicySolver::vertex(DominatingPolicySolver::VertexId id) {
    return VertexOwner::object(id);
}

DominatingPolicySolver::VertexId DominatingPolicySolver::base_vertex() const {
    return VertexOwner::ids().back();
}

double DominatingPolicySolver::rho_for_dominating_vertex(UncertaintyVariable::Index uvar) const {
    return dominating_vertex(uvar)->vertex_value(uvar) - base_vertex()->vertex_value(uvar);
}

double DominatingPolicySolver::approximation_factor() const {
    auto const m = double(model().num_uvars());
    auto const set_type = model().uncertainty_set().special_type();
    switch (mode()) {
        case Mode::THOMAE:
        case Mode::THOMAE_SHIFT: {
            switch (set_type) {
                case UncertaintySet::SpecialSetType::BALL: {
                    return std::pow((std::pow(m, 1. / 2) + 1) / 2, 1. / 2);
                }
                case UncertaintySet::SpecialSetType::BUDGET: {
                    double const b = model().uncertainty_set().budget();
                    helpers::exception_check(b < m, "Budget must be smaller than number of variables!");
                    return b * (m - 1) / (m + b * (b - 2));
                }
                default:
                    helpers::exception_throw("No beta, approximation factor implemented for " +
                                             UncertaintySet::to_string(set_type));
                    return std::nan("0");
            }
        }
        default:
            helpers::exception_throw("Not implemented or not allowed yet!");
            return std::nan("0");
    }
}

SOCVariable::Reference
DominatingPolicySolver::dominating_variable(DecisionVariable::Index const& var, VertexId vertex) const {
    return _dominating_variables[var.raw_id()][vertex.raw_id()];
}

SOCVariable::Reference DominatingPolicySolver::uvar_vertex_shift(UncertaintyVariable::Index uvar) const {
    helpers::exception_check(allow_shift(), "Trying to use uvar shifts when they are disabled!");
    return _vert_shift.at(uvar.raw_id());
}

double DominatingPolicySolver::vertex_value(UncertaintyVariable::Index const& var, VertexId vertex) const {
    return DominatingPolicySolver::vertex(vertex).vertex_value(var);
}

double DominatingPolicySolver::vertex_solution_value(UncertaintyVariable::Index const& var, VertexId vertex) const {
    if (allow_shift()) {
        return uvar_vertex_shift(var)->solution() * (var->ub() - vertex_value(var, vertex)) + vertex_value(var, vertex);
    }
    return vertex_value(var, vertex);
}

AffineExpression<SOCVariable::Reference>
DominatingPolicySolver::vertex_expression(UncertaintyVariable::Index const& var, VertexId vertex) const {
    if (allow_shift()) {
        return uvar_vertex_shift(var) * (var->ub() - vertex_value(var, vertex)) + vertex_value(var, vertex);
    }
    return AffineExpression<SOCVariable::Reference>(vertex_value(var, vertex));
}

void DominatingPolicySolver::build_vertices() {
    helpers::exception_check(num_vertices() == 0, "Vertices mustn't be filled before build!");
    _dominating_vertex_for_uvar = std::vector<VertexId>(model().num_uvars(), VertexId::NO_INDEX);
    auto const set_type = (mode() == Mode::BOX) ?
                          UncertaintySet::SpecialSetType::BOX :
                          model().uncertainty_set().special_type();
    switch (set_type) {
        case UncertaintySet::SpecialSetType::BOX: {
            std::vector<double> v(model().num_uvars());
            for (auto const& var: model().uncertainty_variables()) {
                v[var.id().raw_id()] = var.ub();
            }
            add_vertex(v, model().uncertainty_variable_ids());
            return;
        }
        case UncertaintySet::SpecialSetType::BALL:
        case UncertaintySet::SpecialSetType::BUDGET: {
            switch (mode()) {
                case Mode::BENTAL: {
                    auto const [beta, gamma] = get_beta_gamma_small_approach(model().uncertainty_set().special_type());
                    for (auto const& uvar: model().uncertainty_variables()) {
                        helpers::exception_check(uvar.is_0_1_var(),
                                                 "Dominating Model only defined for uncertainty in[0,1]^d");
                        std::vector<double> v(model().num_uvars(), 0);
                        v[uvar.id().raw_id()] += beta;
                        auto const v_id = add_vertex(v, {uvar.id()});
                    }
                    add_vertex(std::vector<double>(model().num_uvars(), beta * gamma), {});
                    return;
                }
                case Mode::THOMAE:
                case Mode::THOMAE_SHIFT: {
                    auto const [beta, gamma] = get_beta_gamma(model().uncertainty_set().special_type());
                    for (auto const& uvar: model().uncertainty_variables()) {
                        helpers::exception_check(uvar.is_0_1_var(),
                                                 "Dominating Model only defined for uncertainty in[0,1]^d");
                        std::vector<double> v(model().num_uvars(), beta * gamma / 2);
                        v[uvar.id().raw_id()] += beta / 2;
                        auto const v_id = add_vertex(v, {uvar.id()});
                    }
                    add_vertex(std::vector<double>(model().num_uvars(), beta * gamma / 2), {});
                    return;
                }
                default:
                    helpers::exception_throw("Case should not exist!");
            }
            default: {
                helpers::exception_throw("Not implemented yet!");
            }
        }
    }
}

DominatingPolicySolver::VertexId DominatingPolicySolver::add_vertex(std::vector<double> const& vertex,
                                                                    std::vector<UncertaintyVariable::Index> const& dominated_variables) {
    auto const v_id = VertexOwner::base_add_object(*this, vertex, dominated_variables);
    for (auto const uvar: dominated_variables) {
        _dominating_vertex_for_uvar.at(uvar.raw_id()) = v_id;
    }
    return v_id;
}

std::pair<double, double> DominatingPolicySolver::get_beta_gamma(UncertaintySet::SpecialSetType set_type) {
    auto const m = double(model().num_uvars());
    switch (set_type) {
        case UncertaintySet::SpecialSetType::BALL: {
            return {std::pow(m, 1. / 4),
                    1. / std::pow(m, 1. / 2)};
        }
        case UncertaintySet::SpecialSetType::BUDGET: {
            double const b = model().uncertainty_set().budget();
            helpers::exception_check(b < m, "Budget must be smaller than number of variables!");
            return {2 * b * (m - b) / (m + b * (b - 2)),
                    (b - 1) / (m - b)};
        }
        default:
            helpers::exception_check(false,
                                     "No beta, gamma computation allowed for " +
                                     UncertaintySet::to_string(set_type));
            return {};
    }
}

std::pair<double, double>
DominatingPolicySolver::get_beta_gamma_small_approach(UncertaintySet::SpecialSetType set_type) {
    auto const m = double(model().num_uvars());
    switch (set_type) {
        case UncertaintySet::SpecialSetType::BALL: {
            helpers::exception_check(model().uncertainty_set().budget() == 1,
                                     "Ball domination expects a ball of radius 1!"
            );
            return {std::pow(m, 1. / 4),
                    1. / std::pow(m, 1. / 2)};
        }
        case UncertaintySet::SpecialSetType::BUDGET: {
            double const b = model().uncertainty_set().budget();
            helpers::exception_check(b <= m, "Budget must be smaller than number of variables!");
            return {std::min(b, m / b),
                    b / m};
        }
        default:
            helpers::exception_check(false,
                                     "No beta, gamma computation allowed for " +
                                     UncertaintySet::to_string(set_type));
            return {};
    }
}

void DominatingPolicySolver::build_variables() {
    for (auto const& dvar: model().decision_variables()) {
        _dominating_variables.emplace_back(
                soc_model().add_variables(num_vertices(), dvar.name(), dvar.lb(), dvar.ub()));
    }
    if (allow_shift()) {
        for (auto const& uvar: model().uncertainty_variables()) {
            _vert_shift.emplace_back(soc_model().add_variable(uvar.name() + "_shift", 0, 1));
        }
    }
}

void DominatingPolicySolver::build_fix_period_constraints() {
    for (auto const v: VertexOwner::ids()) {
        for (auto const& dvar: model().decision_variables()) {
            if (std::all_of(dvar.dependencies().begin(),
                            dvar.dependencies().end(),
                            [&](VariableReference<UncertaintyVariable> const& uvar) {
                                return vertex_value(uvar, v) == vertex_value(uvar, base_vertex());
                            })) {
                soc_model().add_constraint(
                        dominating_variable(dvar.id(), v) == dominating_variable(dvar.id(), base_vertex()),
                        "Fix " + dvar.name() + " on vertex " + std::to_string(v.raw_id()));
            }
        }
    }
}

void DominatingPolicySolver::build_constraints() {
    for (auto const v: VertexOwner::ids()) {
        for (auto const& constr: model().constraints()) {
            auto const& old_affine_expression = constr.expression();
            soc_model().add_constraint({constr.sense(), ro_affine_to_affine(constr.expression(), v),
                                        constr.name() + std::to_string(v.raw_id())});
        }
    }
}

void DominatingPolicySolver::build_objective() {
    auto const& obj = soc_model().add_variable("Objective", 0, NO_VARIABLE_UB);
    AffineExpression<SOCVariable::Reference> all_vertex_objectives_sum;
    for (auto const v: VertexOwner::ids()) {
        auto obj_vertex = ro_affine_to_affine(model().objective().expression(), v);
        all_vertex_objectives_sum += obj_vertex;
        soc_model().add_constraint(obj_vertex <= obj, "Dominate Vertex " + std::to_string(v.raw_id()));
    }
    soc_model().add_objective({ObjectiveSense::MIN, AffineExpression<SOCVariable::Reference>(obj)});
    if (_reoptimize) {
        soc_model().add_objective({ObjectiveSense::MIN, all_vertex_objectives_sum});
    }
}

std::size_t DominatingPolicySolver::num_vertices() const {
    return VertexOwner::num_objects();
}

AffineExpression<SOCVariable::Reference>
DominatingPolicySolver::ro_affine_to_affine(RoAffineExpression const& old, VertexId v) const {
    AffineExpression<SOCVariable::Reference> new_affine_expression(old.constant());
    for (auto const& scaled_dvar: old.decisions().scaled_variables()) {
        new_affine_expression += ScaledVariable<SOCVariable::Reference>{scaled_dvar.scale(),
                                                                        dominating_variable(
                                                                                scaled_dvar.variable(), v)};
    }
    for (auto const& scaled_uvar: old.uncertainties().scaled_variables()) {
        new_affine_expression += scaled_uvar.scale() * vertex_expression(scaled_uvar.variable(), v);
    }
    return new_affine_expression;
}

std::string DominatingPolicySolver::first_failed_assumption() const {
    if (model().objective().sense() != ObjectiveSense::MIN)
        return "Minimization Problem";
    if (not model().objective().expression().uncertainties().empty())
        return "No Uncertain Objective";
    if (not model().objective().expression().uncertainty_decisions().empty())
        return "Dominating Models can't handle uncertainty scaled decisions!";

    for (auto const& constr: model().constraints()) {
        if (not constr.expression().uncertainty_decisions().empty()) {
            return "Dominating Models can't handle uncertainty scaled decisions!";
        }
        auto const constant = constr.expression().constant();
        if (not((constant <= 0 and constr.sense() == ConstraintSense::GEQ)
                or (constant >= 0 and constr.sense() == ConstraintSense::LEQ)
                or (constant == 0 and constr.sense() == ConstraintSense::EQ)))
            return "b has to be non negative!";

        for (auto const& scaled_uvar: constr.expression().uncertainties().scaled_variables()) {
            if (not((scaled_uvar.scale() <= 0 and constr.sense() == ConstraintSense::GEQ)
                    or (scaled_uvar.scale() >= 0 and constr.sense() == ConstraintSense::LEQ)))
                return "B has to be non negative!";
        };
    }

    if (std::any_of(model().uncertainty_variables().begin(),
                    model().uncertainty_variables().end(),
                    [](auto const& uvar) { return uvar.lb() < 0; }))
        return "h has to be non negative";

    if (std::any_of(model().uncertainty_variables().begin(),
                    model().uncertainty_variables().end(),
                    [](auto const& uvar) { return uvar.ub() == NO_VARIABLE_UB; }))
        return "h has to have upper bounds";

    return "";
}

Vertex DominatingPolicySolver::merge_vertices(std::vector<VertexId> vertices) const {
    std::vector<double> new_vertex(model().num_uvars(), 0);
    std::vector<UncertaintyVariable::Index> dominated_uncertainties;
    for (auto const v: vertices) {
        for (auto const& uvar: model().uncertainty_variables()) {
            new_vertex[uvar.id().raw_id()] = std::max(new_vertex[uvar.id().raw_id()], v->vertex_value(uvar.id()));
        }
        dominated_uncertainties.insert(dominated_uncertainties.end(),
                                       v->dominated_uncertainties().begin(),
                                       v->dominated_uncertainties().end());
    }
    return Vertex(VertexId::NO_INDEX, *this, new_vertex, dominated_uncertainties);
}

void DominatingPolicySolver::override_vertices(std::vector<Vertex> const& new_vertices) {
    VertexOwner::clear();
    for (auto const& v: new_vertices) {
        add_vertex(v.vertex_vector(), v.dominated_uncertainties());
    }
    helpers::exception_check(
            all_uncertainties_dominated_exactly_once(),
            "All uncertainty variables have to be dominated by exactly one vertex for a valid vertex set.");
}

bool DominatingPolicySolver::all_uncertainties_dominated_exactly_once() {
    std::vector<bool> dominated(model().num_uvars(), false);
    for (auto const& v: VertexOwner::objects()) {
        for (auto const& uvar: v.dominated_uncertainties()) {
            if (dominated.at(uvar.raw_id())) {
                return false;
            }
            dominated.at(uvar.raw_id()) = true;
        }
    }
    return std::all_of(dominated.begin(), dominated.end(), [](auto b) { return b; });
}

void DominatingPolicySolver::merge_vertices(std::vector<std::vector<VertexId>> const& vertices_to_merge) {
    std::vector<Vertex> new_vertices;
    std::vector<bool> vertex_merged(num_vertices(), false);
    for (auto const& to_merge: vertices_to_merge) {
        new_vertices.emplace_back(merge_vertices(to_merge));
        for (auto const& v: to_merge) {
            helpers::exception_check(not vertex_merged.at(v.raw_id()), "Each vertex may only be merged once!");
            helpers::exception_check(v.raw_id() != base_vertex().raw_id(), "Do not merge the base vertex!");
            vertex_merged.at(v.raw_id()) = true;
        }
    }
    for (auto const& v: VertexOwner::objects()) {
        if (not vertex_merged.at(v.id().raw_id())) {
            new_vertices.emplace_back(v);
        }
    }
    override_vertices(new_vertices);
}

DominatingPolicySolver::VertexId DominatingPolicySolver::dominating_vertex(UncertaintyVariable::Index uvar) const {
    return _dominating_vertex_for_uvar.at(uvar.raw_id());
}

void DominatingPolicySolver::set_reoptimize(bool reoptimize) {
    _reoptimize = reoptimize;
}

SOCModel& DominatingPolicySolver::soc_model() {
    return _soc_model;
}

SOCModel const& DominatingPolicySolver::soc_model() const {
    return _soc_model;
}

solvers::SOCSolverBase& DominatingPolicySolver::soc_solver() {
    return *_soc_solver;
}

DominatingPolicySolver::Mode const& DominatingPolicySolver::mode() const {
    return _mode;
}

bool DominatingPolicySolver::allow_shift() const {
    return mode() == Mode::THOMAE_SHIFT;
}

}
