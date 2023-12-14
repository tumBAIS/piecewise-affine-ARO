//
// Created by simon on 14/06/2021.
//

#ifndef ROBUSTOPTIMIZATION_CONSTRAINT_H
#define ROBUSTOPTIMIZATION_CONSTRAINT_H

#include "RawConstraint.h"

namespace robust_model {

template<class E>
class Constraint : public RawConstraint<E> {
public:
    using NativeRawConstraint = RawConstraint<E>;

public:
    Constraint(NativeRawConstraint raw_constraint, std::string const& name);

    Constraint(std::string const& name, ConstraintSense sense, E expression);

    template<ConvertibleTo<E> T>
    Constraint(RawConstraint<T> raw_constraint, std::string const& name);

    std::string const& name() const;

    std::string to_string() const override;

    double dual_lb() const;

    double dual_ub() const;

private:
    std::string _name;
};

}

#include "Constraint.tplt"

#endif //ROBUSTOPTIMIZATION_CONSTRAINT_H
