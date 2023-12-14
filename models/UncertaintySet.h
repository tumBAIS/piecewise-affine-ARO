//
// Created by simon on 12.05.21.
//

#ifndef ROBUSTOPTIMIZATION_UNCERTAINTYSET_H
#define ROBUSTOPTIMIZATION_UNCERTAINTYSET_H

#include "basic_model_objects/UncertaintyVariable.h"
#include "basic_model_objects/SOCConstraint.h"
#include <optional>
#include <memory>

namespace robust_model {

class SOCModel;

struct UncertaintyRealization {
public:
    explicit UncertaintyRealization(std::vector<double> values);

    std::vector<double> const& values() const;

    double value(UncertaintyVariable::Index const& id) const;

private:
    std::vector<double> _values;
};

class UncertaintySet;

class UncertaintySetConstraintsSet : public helpers::IndexedObject<UncertaintySetConstraintsSet> {
    using Constraint = SOCConstraint<UncertaintyVariable>;
public:
    UncertaintySetConstraintsSet(Index id, UncertaintySet const& uncertainty_set);

    void add_uncertainty_constraint(Constraint const& constraint);

    std::vector<Constraint> const& constraints() const;

    bool is_box() const;

    bool is_empty() const;

    bool is_norm_ball() const;

    double budget() const;

    void clear();

private:
    UncertaintySet const& _uncertainty_set;
    std::vector<Constraint> _uncertainty_constraints;
};

class ROModel;

class UncertaintySet :
        public helpers::IndexedObjectOwner<UncertaintyVariable>,
        private helpers::IndexedObjectOwner<UncertaintySetConstraintsSet> {
    friend UncertaintyRealization;
public:
    using UncertaintyReference = VariableReference<UncertaintyVariable>;
    using Constraint = SOCConstraint<UncertaintyVariable>;

public:
    enum class SpecialSetType {
        BOX,
        BALL,
        BUDGET,
        OTHER
    };

    static std::string to_string(SpecialSetType set_type);

public:
    explicit UncertaintySet(ROModel const& model);

    UncertaintyReference add_variable(std::string const& name,
                                      std::optional<period_id> p = std::optional<period_id>(),
                                      double lb = NO_VARIABLE_LB,
                                      double ub = NO_VARIABLE_UB);

    std::vector<UncertaintyVariable> const& variables() const;

    std::vector<double> lower_bounds() const;

    std::vector<double> upper_bounds() const;

    size_t num_variables() const;

    std::vector<UncertaintyVariable::Index> const& indices() const;

    ROModel const& model() const;

    SpecialSetType special_type() const;

    void add_special_type_constraint(SpecialSetType type, double budget=1);

    // in case of multiple constraint sets this will add the constraint to last set!
    void add_uncertainty_constraint(Constraint const& constraint);

    void add_uncertainty_constraint(Constraint const& constraint, UncertaintySetConstraintsSet::Index constraint_set);

    std::vector<Constraint> const& uncertainty_constraints() const;

    std::vector<Constraint> const& uncertainty_constraints(UncertaintySetConstraintsSet::Index constraint_set) const;

    UncertaintySetConstraintsSet::Index add_constraint_set();

    std::vector<UncertaintySetConstraintsSet::Index> const& constraint_sets() const;

    // Warning this does not invalidate old smart indices!
    // They might simply refer to a new object or exceed the vector length!
    void clear_constraints();

    std::string full_string() const;

    std::unique_ptr<SOCModel> to_soc_model() const;

    double budget() const;

    bool non_negative() const;

    bool all_0_1_variables() const;

    bool all_bounded_variables() const;

    // note that free sets are not checked for rotational invariance -> they return false
    bool rotational_invariant() const;

    // note that free sets are not checked for symmetry -> they return false
    bool symmetric() const;

    double max_one_norm_k_active(size_t k) const;

    std::vector<std::vector<double>> generate_uncertainty(size_t num_realizations) const;

private:
    std::vector<std::vector<double>> generate_uncertainty_ball(size_t num_realizations) const;

    std::vector<std::vector<double>> generate_uncertainty_budgeted(size_t num_realizations) const;

private:
    ROModel const& _model;
};

}

#endif //ROBUSTOPTIMIZATION_UNCERTAINTYSET_H
