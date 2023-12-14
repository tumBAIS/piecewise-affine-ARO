//
// Created by simon on 11.05.21.
//

#ifndef ROBUSTOPTIMIZATION_SOCCONSTRAINT_H
#define ROBUSTOPTIMIZATION_SOCCONSTRAINT_H

#include "SOCExpression.h"
#include "Constraint.h"
#include "../../helpers/helpers.h"

namespace robust_model {

template<class V>
class SOCConstraint : public Constraint<SOCExpression<V>> {
public:
    using BaseConstraint = Constraint<SOCExpression<V>>;
    using RawConstraint = typename BaseConstraint::RawConstraint;

public:

    SOCConstraint(RawConstraint const& raw_constraint, std::string const& name);

    template<ConvertibleTo<typename SOCConstraint<V>::RawConstraint> RC>
    SOCConstraint(RC const& raw_constraint, std::string const& name);

    template<ConvertibleTo<SOCExpression<V>> E>
    SOCConstraint(ConstraintSense sense, E const& expression, std::string const& name);

    template<class W>
    SOCConstraint<W> translate_to_other(std::vector<VariableReference<W> > const& new_vars) const;

    template<class W, AffineAddable<typename W::Reference> T>
    SOCConstraint<W> substitute(std::vector<T> const& substitutions) const;

    SOCExpression<V> const& soc_expression() const;

    bool is_variable_bound() const;

    bool has_dual_value() const;

    double dual_value() const;

    void set_dual_value(double value);
    void invalidate_dual_value();

    template<class S>
    bool constraint_satisfied(S const& realization) const;

private:
    std::optional<double> _dual_value;

};

}

#include "SOCConstraint.tplt"

#endif //ROBUSTOPTIMIZATION_SOCCONSTRAINT_H
