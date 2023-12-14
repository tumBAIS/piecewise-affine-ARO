//
// Created by thomae on 22/06/2023.
//

#ifndef ROBUSTOPTIMIZATION_BGPIECEWISELINEARPOLICYSOLVER_H
#define ROBUSTOPTIMIZATION_BGPIECEWISELINEARPOLICYSOLVER_H

#include "../../models/ROModel.h"
#include "../../models/SOCModel.h"
#include "AROPolicySolverBase.h"
#include "../soc_solvers/GurobiSOCSolver.h"


namespace robust_model {

class BGPiecewiseLinearPolicy {
public:
    BGPiecewiseLinearPolicy(DecisionVariable::Index const& decision_variable,
                            std::vector<std::vector<double>> x_upper,
                            std::vector<std::vector<double>> x_lower,
                            std::vector<double> constants_upper,
                            std::vector<double> constants_lower);

    double value(std::vector<double> const& uncertainty_realization) const;

    std::string to_string() const;

private:
    double find_max(
            std::vector<double> const& uncertainty_realization,
            std::vector<std::vector<double>> const& scales,
            std::vector<double> const& constants) const;

    std::string max_string(
            std::vector<std::vector<double>> const& scales,
            std::vector<double> const& constants) const;

private:
    DecisionVariable::Index const _decision_variable;
    std::vector<std::vector<double>> const _x_upper;
    std::vector<std::vector<double>> const _x_lower;
    std::vector<double> const _constants_upper;
    std::vector<double> const _constants_lower;
};

class BGPiecewiseLinearPolicySolver : public solvers::AROPolicySolverBase{
    enum class PiecewisePart {
        UPPER,
        LOWER
    };
    static std::string to_string(PiecewisePart part){
        switch (part) {
            case PiecewisePart::UPPER:
                return "Upper";
            case PiecewisePart::LOWER:
                return "Lower";
        }
    }

    enum class ProblemType {
        MASTER,
        SUBPROBLEM
    };

    template<ProblemType V>
    struct ProblemTypeDecisionMap {
        using Type = void;
    };

    template<ProblemType V>
    struct ProblemTypeUncertaintyMap {
        using Type = void;
    };

public:
    BGPiecewiseLinearPolicySolver(ROModel const& model, size_t num_pieces);

    SolutionRealization specific_solution(std::vector<double> const& uncertainty_realization) const final;

    std::string const& name() const;

    std::string current_solution_string() const;

private:
    void build_implementation() final;

    void solve_implementation() final;

    template<ProblemType T>
    SOCModel& soc_model();

    bool sufficiently_feasible() const;

    void add_maximally_violated_scenarios();

    std::tuple<size_t, std::vector<size_t>> find_violations(
            std::vector<std::vector<double>> const& uncertainty_realizations) const;

    bool check_variable_violations(
            size_t scenario,
            SolutionRealization const& solution_realization,
            std::vector<double>& var_bound_violation_amount,
            std::vector<size_t>& var_bound_violation_scenario) const;

    bool check_constraint_violations(
            size_t scenario,
            SolutionRealization const& solution_realization,
            std::vector<double>& constraint_violation_amount,
            std::vector<size_t>& constraint_violation_scenario) const;

    bool check_objective_violations(
            size_t scenario,
            SolutionRealization const& solution_realization,
            double& objective_violation_amount,
            size_t& objective_violation_scenario) const;


    void add_base_model_variables();

    template<ProblemType T>
    typename ProblemTypeDecisionMap<T>::Type
    piecewise_variable(DecisionVariable::Index const& dvar, size_t piece, size_t dependency_pos,
                       PiecewisePart part) const;

    template<ProblemType T>
    typename ProblemTypeDecisionMap<T>::Type
    piecewise_constant(DecisionVariable::Index const& dvar, size_t piece, PiecewisePart part) const;

    void add_scenario(std::vector<double> const& uncertainty_realization);

    template<ProblemType T>
    std::vector<SOCVariable::Reference> add_piecewise_function_encoding(
            std::vector<typename ProblemTypeUncertaintyMap<T>::Type> const& uncertainty);

    template<ProblemType T>
    SOCVariable::Reference add_maximum_variable_with_constraints(
            DecisionVariable::Index const& dvar,
            std::vector<typename ProblemTypeUncertaintyMap<T>::Type> const& uncertainty_realization,
            PiecewisePart part, std::string const& var_name_addendum = "");

    void add_scenario_constraints(std::vector<double> const& uncertainty_realization,
                                  std::vector<SOCVariable::Reference> const& scenario_var);

    void add_scenario_objective(std::vector<double> const& uncertainty_realization,
                                std::vector<SOCVariable::Reference> const& scenario_var);

    void extract_current_solution();

    size_t calculate_verification_sample_size() const;

    size_t calculate_training_sample_size() const;

    void build_submodel();

    // functions return true, if a constraint was added and false if not
    bool add_maximally_violating_uncertainty_realizations_with_submodel();

    bool add_violating_uncertainty_for_maximizing_expression(RoAffineExpression const& expression);

    solvers::SOCSolverBase & soc_solver();
    solvers::SOCSolverBase const& soc_solver() const;

    solvers::GurobiSOCSolver & submodel_soc_solver();

private:
    size_t const _num_pieces;

    double const _delta = 0.005;
    double const _beta = 1e-4;
    double const _beta_hat = 1e-4;
    double const _eps_hat = 0.005;
    size_t const _verification_sample_size;
    size_t const _training_sample_size;

    // These are the upper and lower bounds for the factors of the piecewise affine functions
    // They are needed for numeric stability
    double const _var_bounds = 10;

    std::vector<std::vector<double>> _training_samples;

    double const _tolerance = 1e-2;

    SOCModel _solver_model;
    std::unique_ptr<solvers::GurobiSOCSolver> _soc_solver;

    std::unique_ptr<SOCModel> _submodel;
    std::unique_ptr<solvers::GurobiSOCSolver> _submodel_soc_solver;

    std::unique_ptr<std::vector<SOCVariable::Reference>> _submodel_uvars;
    std::unique_ptr<std::vector<SOCVariable::Reference>> _submodel_dvars;

    // format of piecewise affine policy variables:
    // variables x pieces x dependent_uncertainties
    std::vector<std::vector<std::vector<SOCVariable::Reference>>> _upper_variables;
    std::vector<std::vector<std::vector<SOCVariable::Reference>>> _lower_variables;
    // format of piecewise affine policy constants variables:
    // variables x pieces
    std::vector<std::vector<SOCVariable::Reference>> _upper_variables_constant;
    std::vector<std::vector<SOCVariable::Reference>> _lower_variables_constant;

    std::vector<BGPiecewiseLinearPolicy> _current_solution;
};

}

#include "BGPiecewiseLinearPolicySolver.tplt"

#endif //ROBUSTOPTIMIZATION_BGPIECEWISELINEARPOLICYSOLVER_H
