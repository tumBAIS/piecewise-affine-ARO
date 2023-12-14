//
// Created by thomae on 22/06/2023.
//

#include "BGPiecewiseLinearPolicySolver.h"

#include <utility>
#include <set>
#include <cmath>

namespace robust_model {

BGPiecewiseLinearPolicy::BGPiecewiseLinearPolicy(
        DecisionVariable::Index const& decision_variable,
        std::vector<std::vector<double>> x_upper,
        std::vector<std::vector<double>> x_lower,
        std::vector<double> constants_upper,
        std::vector<double> constants_lower) :
        _decision_variable(decision_variable),
        _x_upper(std::move(x_upper)),
        _x_lower(std::move(x_lower)),
        _constants_upper(std::move(constants_upper)),
        _constants_lower(std::move(constants_lower)) {}

double BGPiecewiseLinearPolicy::value(std::vector<double> const& uncertainty_realization) const {
    double const upper_max = find_max(uncertainty_realization, _x_upper, _constants_upper);
    double const lower_max = find_max(uncertainty_realization, _x_lower, _constants_lower);
    return upper_max - lower_max;
}

double BGPiecewiseLinearPolicy::find_max(
        std::vector<double> const& uncertainty_realization,
        std::vector<std::vector<double>> const& scales,
        std::vector<double> const& constants) const {
    double maximum = -std::numeric_limits<double>::max();
    for (size_t piece = 0; piece < scales.size(); ++piece) {
        double val = constants.at(piece);
        for (size_t dependency = 0; dependency < _decision_variable->dependencies().size(); ++dependency) {
            val += scales.at(piece).at(dependency) *
                   uncertainty_realization.at(_decision_variable->dependencies().at(dependency).raw_id());
        }
        maximum = std::max(maximum, val);
    }
    return maximum;
}

std::string BGPiecewiseLinearPolicy::to_string() const {
    return _decision_variable->name() + "=" +
           max_string(_x_upper, _constants_upper) + " - \n" +
           max_string(_x_lower, _constants_lower);
}

std::string BGPiecewiseLinearPolicy::max_string(std::vector<std::vector<double>> const& scales,
                                                std::vector<double> const& constants) const {
    std::string s = "MAX(";
    for (size_t piece = 0; piece < constants.size(); ++piece) {
        s += "\n";
        for (size_t dependency = 0; dependency < scales.at(piece).size(); ++dependency) {
            s += std::to_string(scales.at(piece).at(dependency)) + "*" +
                 _decision_variable->dependencies().at(dependency)->name() + " + ";
        }
        s += std::to_string(constants.at(piece));
    }
    s += ")";
    return s;
}

BGPiecewiseLinearPolicySolver::BGPiecewiseLinearPolicySolver(ROModel const& model, size_t const num_pieces) :
        AROPolicySolverBase(model), _num_pieces(num_pieces),
        _verification_sample_size(calculate_verification_sample_size()),
        _training_sample_size(calculate_training_sample_size()),
        _solver_model("BGPAPModel" + model.name()) {}


SolutionRealization
BGPiecewiseLinearPolicySolver::specific_solution(std::vector<double> const& uncertainty_realization) const {
    std::vector<double> solutions(model().num_dvars());
    for (auto const& dvar: model().decision_variable_ids()) {
        solutions.at(dvar.raw_id()) = _current_solution.at(dvar.raw_id()).value(uncertainty_realization);
    }
    return {model(), uncertainty_realization, std::move(solutions)};
}

void BGPiecewiseLinearPolicySolver::build_implementation() {
    add_base_model_variables();
    auto const obj_variable = _solver_model.add_variable("Obj");
    _solver_model.add_objective({ObjectiveSense::MIN, AffineExpression<SOCVariable::Reference>(obj_variable)});
    add_scenario(model().uncertainty_set().generate_uncertainty(1).at(0));
//    _training_samples = model().uncertainty_set().generate_uncertainty(_training_sample_size);
    _soc_solver = std::make_unique<solvers::GurobiSOCSolver>(soc_model<ProblemType::MASTER>());
    set_parameters_to_other(soc_solver());
}

std::string const& BGPiecewiseLinearPolicySolver::name() const {
    return _solver_model.name();
}

//void BGPiecewiseLinearPolicySolver::solve_with_gurobi() {
//    auto start = std::chrono::high_resolution_clock::now();
//    auto const current_runtime = [&]() {
//        return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
//    };
//    double model_time = 0;
//    _solver_model.set_maximal_runtime(_maximal_runtime - current_runtime());
//    _solver_model.solve_with_gurobi();
//    model_time += _solver_model.solution_time();
//    extract_current_solution();
//    while (not sufficiently_feasible()) {
//        add_maximally_violated_scenarios();
//        _solver_model.set_maximal_runtime(_maximal_runtime - current_runtime());
//        _solver_model.solve_with_gurobi();
//        model_time += _solver_model.solution_time();
//        if (current_runtime() > _maximal_runtime) {
//            break;
//        }
//        extract_current_solution();
//    }
//    _solution_time = (current_runtime() <= _maximal_runtime) ? current_runtime() : -1.;
//}


void BGPiecewiseLinearPolicySolver::solve_implementation() {
    set_parameters_to_other(soc_solver());
    auto start = std::chrono::high_resolution_clock::now();
    auto const current_runtime = [&]() {
        return std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start).count();
    };
    double model_time = 0;

    auto const solve_master_problem = [&]() {
        if(has_runtime_limit())
            soc_solver().set_runtime_limit(runtime_limit() - current_runtime());
        soc_solver().solve();
        if (not soc_solver().has_solution()) {
            return false;
        }
        model_time += soc_solver().runtime();
        extract_current_solution();
        return true;
    };

    if (not solve_master_problem())
        return;

    while (add_maximally_violating_uncertainty_realizations_with_submodel()) {
        if (has_runtime_limit() and current_runtime() > runtime_limit() - 1.) {
            break;
        }
        if (not solve_master_problem())
            return;
    }
    double const total_runtime = current_runtime();
    if (not has_runtime_limit() or total_runtime <= runtime_limit() - 1.) {
        set_status(solvers::SolverBase::Status::OPTIMAL);
        set_objective_value(soc_solver().objective_value());
        set_runtime(total_runtime);
    } else {
        set_status(solvers::SolverBase::Status::TIME_LIMIT);
    }
}

bool BGPiecewiseLinearPolicySolver::sufficiently_feasible() const {
    auto const [num_violations, violated_scenarios] = find_violations(
            model().uncertainty_set().generate_uncertainty(_verification_sample_size)
    );
    return double(num_violations) / double(_verification_sample_size) <= _delta;
}

void BGPiecewiseLinearPolicySolver::add_maximally_violated_scenarios() {
    auto const [num_violations, violated_scenarios] = find_violations(
            _training_samples
    );
    for (auto const& scenario_id: violated_scenarios) {
        add_scenario(_training_samples.at(scenario_id));
    }
}

#define TIMEPROCESS(PROCESS, TIMEVAR) \
        auto start_##TIMEVAR = std::chrono::high_resolution_clock::now(); \
        PROCESS;               \
        TIMEVAR += std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - start_##TIMEVAR).count();

std::tuple<size_t, std::vector<size_t>>
BGPiecewiseLinearPolicySolver::find_violations(std::vector<std::vector<double>> const& uncertainty_realizations) const {
    size_t num_violated_scenarios = 0;
    std::vector<double> constraint_violation_amount(model().constraints().size(), 0.);
    std::vector<size_t> constraint_violation_scenario(model().constraints().size());
    std::vector<double> var_bound_violation_amount(model().num_dvars(), 0.);
    std::vector<size_t> var_bound_violation_scenario(model().num_dvars());
    double objective_violation_amount(0.);
    size_t objective_violation_scenario;
    double solution_creation = 0;
    double solution_evaluation = 0;

    for (size_t scenario = 0; scenario < uncertainty_realizations.size(); ++scenario) {
        auto const& realization = uncertainty_realizations.at(scenario);
        TIMEPROCESS(auto const solution_realization = specific_solution(realization), solution_creation)
        TIMEPROCESS(
                bool const current_scenario_violated =
                        check_variable_violations(scenario, solution_realization,
                                                  var_bound_violation_amount, var_bound_violation_scenario) or
                        check_constraint_violations(scenario, solution_realization,
                                                    constraint_violation_amount, constraint_violation_scenario) or
                        check_objective_violations(scenario, solution_realization,
                objective_violation_amount, objective_violation_scenario), solution_evaluation)

        if (current_scenario_violated) {
            ++num_violated_scenarios;
        }
    }

    //Create unique set of maximally violated scenarios
    std::set<size_t> maximally_violated_scenarios;
    auto const add_to_violated_scenarios = [&](
            std::vector<double> const& violation_amount,
            std::vector<size_t> const& violation_scenario) {
        for (size_t i = 0; i < violation_amount.size(); ++i) {
            if (violation_amount.at(i) > 0) {
                maximally_violated_scenarios.emplace(violation_scenario.at(i));
            }
        }
    };
    add_to_violated_scenarios(constraint_violation_amount, constraint_violation_scenario);
    add_to_violated_scenarios(var_bound_violation_amount, var_bound_violation_scenario);
    if (objective_violation_amount > 0) {
        maximally_violated_scenarios.emplace(objective_violation_scenario);
    }

    return {num_violated_scenarios,
            std::vector<size_t>(maximally_violated_scenarios.begin(),
                                maximally_violated_scenarios.end())};
}

bool BGPiecewiseLinearPolicySolver::check_variable_violations(
        size_t const scenario,
        SolutionRealization const& solution_realization,
        std::vector<double>& var_bound_violation_amount,
        std::vector<size_t>& var_bound_violation_scenario) const {
    bool current_scenario_violated = false;
    for (auto const& dvar: model().decision_variable_ids()) {
        double const bound_violation = std::max(
                dvar->lb() - solution_realization.value(dvar),
                solution_realization.value(dvar) - dvar->ub());
        if (bound_violation > _tolerance) {
            current_scenario_violated = true;
            if (var_bound_violation_amount.at(dvar.raw_id()) < bound_violation) {
                var_bound_violation_amount.at(dvar.raw_id()) = bound_violation;
                var_bound_violation_scenario.at(dvar.raw_id()) = scenario;
            }
        }
    }
    return current_scenario_violated;
}

bool BGPiecewiseLinearPolicySolver::check_constraint_violations(
        size_t const scenario,
        SolutionRealization const& solution_realization,
        std::vector<double>& constraint_violation_amount,
        std::vector<size_t>& constraint_violation_scenario) const {
    bool current_scenario_violated = false;
    for (size_t constraint_id = 0; constraint_id < model().constraints().size(); ++constraint_id) {
        auto const& constraint = model().constraints().at(constraint_id);
        double const lhs = constraint.expression().value(solution_realization);
        double const violation = [&]() {
            switch (constraint.sense()) {
                case ConstraintSense::GEQ:
                    return -lhs;
                case ConstraintSense::LEQ:
                    return lhs;
                case ConstraintSense::EQ:
                    return std::abs(lhs);
            }
        }();
        if (violation > _tolerance) {
            current_scenario_violated = true;
            if (constraint_violation_amount.at(constraint_id) < violation) {
                constraint_violation_amount.at(constraint_id) = violation;
                constraint_violation_scenario.at(constraint_id) = scenario;
            }
        }
    }
    return current_scenario_violated;
}

bool BGPiecewiseLinearPolicySolver::check_objective_violations(
        size_t const scenario,
        SolutionRealization const& solution_realization,
        double& objective_violation_amount,
        size_t& objective_violation_scenario) const {
    double obj_violation =
            ((model().objective().sense() == ObjectiveSense::MIN) ? 1. : -1) *
            (model().objective_value_for_solution(solution_realization) - soc_solver().objective_value());
    if (obj_violation > _tolerance) {
        if (obj_violation > objective_violation_amount) {
            objective_violation_amount = obj_violation;
            objective_violation_scenario = scenario;
        }
        return true;
    }
    return false;
}

void BGPiecewiseLinearPolicySolver::add_base_model_variables() {
    for (auto const& dvar: model().decision_variables()) {
        _upper_variables.emplace_back();
        _lower_variables.emplace_back();
        _upper_variables_constant.emplace_back();
        _lower_variables_constant.emplace_back();
        for (size_t piece = 0; piece < _num_pieces; ++piece) {
            _upper_variables.back().emplace_back();
            _lower_variables.back().emplace_back();
            _upper_variables_constant.back().emplace_back(
                    _solver_model.add_variable(dvar.name() + "_UC_" + std::to_string(piece),
                                               -_var_bounds, _var_bounds));
            _lower_variables_constant.back().emplace_back(
                    _solver_model.add_variable(dvar.name() + "_LC_" + std::to_string(piece),
                                               -_var_bounds, _var_bounds));
            for (auto const& uvar: dvar.dependencies()) {
                _upper_variables.back().back().emplace_back(
                        _solver_model.add_variable(dvar.name() + "_UV_" + std::to_string(piece) + "_" + uvar->name(),
                                                   -_var_bounds, _var_bounds));
                _lower_variables.back().back().emplace_back(
                        _solver_model.add_variable(dvar.name() + "_LV_" + std::to_string(piece) + "_" + uvar->name(),
                                                   -_var_bounds, _var_bounds));
            }
        }
    }
}


void BGPiecewiseLinearPolicySolver::add_scenario(std::vector<double> const& uncertainty_realization) {
    auto const scenario_vars = add_piecewise_function_encoding<ProblemType::MASTER>(uncertainty_realization);
    add_scenario_constraints(uncertainty_realization, scenario_vars);
    add_scenario_objective(uncertainty_realization, scenario_vars);
}

void BGPiecewiseLinearPolicySolver::extract_current_solution() {
    _current_solution.clear();
    for (auto const& dvar: model().decision_variable_ids()) {
        std::vector<std::vector<double>> x_upper(_num_pieces, std::vector<double>(dvar->dependencies().size()));
        std::vector<std::vector<double>> x_lower(_num_pieces, std::vector<double>(dvar->dependencies().size()));
        std::vector<double> constants_upper(_num_pieces);
        std::vector<double> constants_lower(_num_pieces);

        for (size_t piece = 0; piece < _num_pieces; ++piece) {
            for (size_t dependency_pos = 0; dependency_pos < dvar->dependencies().size(); ++dependency_pos) {
                x_upper.at(piece).at(dependency_pos) =
                        piecewise_variable<ProblemType::SUBPROBLEM>(dvar, piece, dependency_pos, PiecewisePart::UPPER);
                x_lower.at(piece).at(dependency_pos) =
                        piecewise_variable<ProblemType::SUBPROBLEM>(dvar, piece, dependency_pos, PiecewisePart::LOWER);
            }
            constants_upper.at(piece) = piecewise_constant<ProblemType::SUBPROBLEM>(dvar, piece, PiecewisePart::UPPER);
            constants_lower.at(piece) = piecewise_constant<ProblemType::SUBPROBLEM>(dvar, piece, PiecewisePart::LOWER);
        }
        _current_solution.emplace_back(dvar, std::move(x_upper), std::move(x_lower),
                                       std::move(constants_upper), std::move(constants_lower));
    }
}

void BGPiecewiseLinearPolicySolver::add_scenario_constraints(std::vector<double> const& uncertainty_realization,
                                                             std::vector<SOCVariable::Reference> const& scenario_var) {
    for (auto const& constraint: model().constraints()) {
        _solver_model.add_constraint(
                {constraint.sense(),
                 constraint.expression().substitute_to_other_affine<typename SOCVariable::Reference>(
                         scenario_var,
                         uncertainty_realization),
                 "ScenarioConstraint" + constraint.name()});
    }
}

void BGPiecewiseLinearPolicySolver::add_scenario_objective(std::vector<double> const& uncertainty_realization,
                                                           std::vector<SOCVariable::Reference> const& scenario_var) {
    _solver_model.add_constraint({
                                         (model().objective().sense() == ObjectiveSense::MIN) ? ConstraintSense::LEQ
                                                                                             : ConstraintSense::GEQ,
                                         model().objective().expression().substitute_to_other_affine<typename SOCVariable::Reference>(
                                                 scenario_var,
                                                 uncertainty_realization)
                                         - _solver_model.objective().expression(),
                                         "ScenarioObjectiveBound"
                                 });
}

size_t BGPiecewiseLinearPolicySolver::calculate_verification_sample_size() const {
    return size_t(std::log10(2. / _beta_hat) / (2. * _eps_hat * _eps_hat)) + 1;
}

size_t BGPiecewiseLinearPolicySolver::calculate_training_sample_size() const {
    double V = 2. *
               double(1 + 2 * _num_pieces * model().num_dvars() * model().num_uvars()) *
               std::log2(8. * std::exp(1) *
                         double(model().constraints().size() + 1) *
                         double(_num_pieces * model().num_dvars()));
    return size_t(4. / _delta * (V * std::log(12. / _delta) + std::log(2. / _beta)) + 1);
}

void BGPiecewiseLinearPolicySolver::build_submodel() {
    _submodel = model().uncertainty_set().to_soc_model();
    _submodel_uvars = std::make_unique<std::vector<SOCVariable::Reference>>();
    for (auto const& uvar: _submodel->variables()) {
        _submodel_uvars->emplace_back(uvar.reference());
    }
    _submodel_dvars = std::make_unique<std::vector<SOCVariable::Reference>>(
            add_piecewise_function_encoding<ProblemType::SUBPROBLEM>(*_submodel_uvars));
    _submodel_soc_solver = std::make_unique<solvers::GurobiSOCSolver>(soc_model<ProblemType::SUBPROBLEM>());
    set_parameters_to_other(submodel_soc_solver());
}

bool BGPiecewiseLinearPolicySolver::add_maximally_violating_uncertainty_realizations_with_submodel() {
    build_submodel();
    bool scenario_added = false;
    for (auto const& constr: model().constraints()) {
        if (constr.sense() == ConstraintSense::LEQ or constr.sense() == ConstraintSense::EQ) {
            scenario_added = scenario_added or
                             add_violating_uncertainty_for_maximizing_expression(constr.expression());
        }
        if (constr.sense() == ConstraintSense::GEQ or constr.sense() == ConstraintSense::EQ) {
            scenario_added = scenario_added or
                             add_violating_uncertainty_for_maximizing_expression(-constr.expression());
        }
    }
    for (auto const& var: model().decision_variables()) {
        if (var.ub() < NO_VARIABLE_UB) {
            scenario_added = scenario_added or
                             add_violating_uncertainty_for_maximizing_expression(
                                     RoAffineExpression(var.reference() - var.ub()));
        }
        if (var.lb() > NO_VARIABLE_LB) {
            scenario_added = scenario_added or
                             add_violating_uncertainty_for_maximizing_expression(
                                     RoAffineExpression(var.lb() - var.reference()));
        }
    }
    scenario_added = scenario_added or
                     add_violating_uncertainty_for_maximizing_expression(
                             ((model().objective().sense() == ObjectiveSense::MIN) ? 1. : -1.) *
                             (
                                     model().objective().expression() - soc_solver().objective_value()
                             )
                     );
    return scenario_added;
}

bool BGPiecewiseLinearPolicySolver::add_violating_uncertainty_for_maximizing_expression(RoAffineExpression const& expression) {
    soc_model<ProblemType::SUBPROBLEM>().clear_and_set_objective(
            {ObjectiveSense::MAX,
             expression.substitute_to_other_affine_no_uncertain_decisions<SOCVariable::Reference>(
                     *_submodel_dvars,
                     *_submodel_uvars
             )});
    submodel_soc_solver().objectives_reset();
    submodel_soc_solver().solve();

    if (submodel_soc_solver().objective_value() > _tolerance) {
        std::vector<double> realization(_submodel_uvars->size());
        for (size_t decision_id = 0; decision_id < _submodel_uvars->size(); ++decision_id) {
            realization.at(decision_id) = _submodel_uvars->at(decision_id)->solution();
        }
        add_scenario(realization);
        return true;
    }
    return false;
}

std::string BGPiecewiseLinearPolicySolver::current_solution_string() const {
    std::string s;
    for (auto const& var_sol: _current_solution) {
        s += var_sol.to_string() + "\n";
    }
    return s;
}

solvers::SOCSolverBase& BGPiecewiseLinearPolicySolver::soc_solver() {
    return *_soc_solver;
}

solvers::SOCSolverBase const& BGPiecewiseLinearPolicySolver::soc_solver() const{
    return *_soc_solver;
}

solvers::GurobiSOCSolver& BGPiecewiseLinearPolicySolver::submodel_soc_solver() {
    return *_submodel_soc_solver;
}

}
