//
// Created by simon on 11.05.21.
//

#ifndef ROBUSTOPTIMIZATION_ROMODEL_H
#define ROBUSTOPTIMIZATION_ROMODEL_H

#include "basic_model_objects/RoAffineExpression.h"
#include "basic_model_objects/ObjectiveBase.h"
#include "UncertaintySet.h"
#include "../helpers/helpers.h"

namespace robust_model {

class ROModel : public helpers::IndexedObjectOwner<DecisionVariable> {
public:
    using RoConstraint = Constraint<RoAffineExpression>;
    using DecisionReference = VariableReference<DecisionVariable>;
    using UncertaintyReference = UncertaintySet::UncertaintyReference;
    using Objective = ObjectiveBase<RoAffineExpression>;

public:
    explicit ROModel(std::string name = "Model");

    std::vector<DecisionReference>
    add_decision_variables(std::size_t n, std::string const& name, std::optional<period_id> p={},
                           double lb = NO_VARIABLE_LB, double ub = NO_VARIABLE_UB);

    std::vector<DecisionReference>
    add_decision_variables_for_each_period(size_t n, std::string const& name, period_id first_period,
                                           double lb = NO_VARIABLE_LB, double ub = NO_VARIABLE_UB);

    DecisionReference
    add_decision_variable(std::string const& name, std::optional<period_id> p={},
                          double lb = NO_VARIABLE_LB, double ub = NO_VARIABLE_UB);

    DecisionReference
    add_decision_variable(std::string const& name, std::vector<UncertaintyVariable::Reference> const& dependencies,
                          double lb = NO_VARIABLE_LB, double ub = NO_VARIABLE_UB);

    std::vector<UncertaintyReference>
    add_uncertainty_variables(std::size_t n, std::string const& name, std::optional<period_id> p={}, double lb = NO_VARIABLE_LB,
                              double ub = NO_VARIABLE_UB);

    std::vector<UncertaintyReference>
    add_uncertainty_variables_for_each_period(size_t n, std::string const& name, period_id first_period,
                                              double lb = NO_VARIABLE_LB, double ub = NO_VARIABLE_UB);

    UncertaintyReference
    add_uncertainty_variable(std::string const& name, std::optional<period_id> p={},
                             double lb = NO_VARIABLE_LB, double ub = NO_VARIABLE_UB);

    void add_constraint(RoConstraint const& constraint);

    template<ConvertibleTo<RoConstraint::RawConstraint> RC>
    void add_constraint(RC const& raw_constraint, std::string const& name) {
        add_constraint({raw_constraint, name});
    }


    void add_uncertainty_constraint(UncertaintySet::Constraint const& constraint);

    template<ConvertibleTo<UncertaintySet::Constraint::RawConstraint> RC>
    void add_uncertainty_constraint(RC const& raw_constraint, std::string const& name) {
        add_uncertainty_constraint({raw_constraint, name});
    }

    UncertaintySetConstraintsSet::Index add_uncertainty_constraint_set();

    void add_uncertainty_constraint(UncertaintySet::Constraint const& constraint,
                                    UncertaintySetConstraintsSet::Index union_set);

    template<ConvertibleTo<UncertaintySet::Constraint::RawConstraint> RC>
    void add_uncertainty_constraint(RC const& raw_constraint, std::string const& name,
                                    UncertaintySetConstraintsSet::Index union_set) {
        add_uncertainty_constraint({raw_constraint, name}, union_set);
    }


    void set_objective(RoAffineExpression const& objective, Objective::Sense sense);

    std::string full_string() const;

    std::vector<DecisionVariable> const& decision_variables() const;

    std::vector<DecisionVariable::Index> const& decision_variable_ids() const;

    std::vector<UncertaintyVariable> const& uncertainty_variables() const;

    std::vector<UncertaintyVariable::Index> const& uncertainty_variable_ids() const;

    std::size_t num_dvars() const;

    std::size_t num_uvars() const;

    UncertaintySet const& uncertainty_set() const;

    void add_special_uncertainty_constraint(UncertaintySet::SpecialSetType type, double budget=1);

    Objective const& objective() const;

    double objective_value_for_solution(SolutionRealization const& solution) const;

    bool solution_feasible(SolutionRealization const& solution, double tolerance = 1e-3) const;

    std::vector<RoConstraint> const& constraints() const;

    std::string const& name() const;

private:
    UncertaintySet _uncertainty_set;
    std::optional<Objective> _objective;
    std::vector<RoConstraint> _constraints;
    std::string const _name;
};

}

#endif //ROBUSTOPTIMIZATION_ROMODEL_H
