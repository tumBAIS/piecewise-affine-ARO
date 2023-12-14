//
// Created by simon on 14/06/2021.
//

#ifndef ROBUSTOPTIMIZATION_RAWCONSTRAINT_H
#define ROBUSTOPTIMIZATION_RAWCONSTRAINT_H

#include <string>

#include "types_and_constants.h"
#include "predeclarations_and_concepts.h"

namespace robust_model {

template<class E>
class RawConstraint {
public:
    RawConstraint(ConstraintSense sense, E const& expression);

    template<ConvertibleTo<E> T>
    RawConstraint(ConstraintSense sense, T const& expression);

    template<ConvertibleTo<E> T>
    explicit RawConstraint(RawConstraint<T> const& raw_constraint);

    virtual std::string to_string() const;

    ConstraintSense sense() const;

    E const& expression() const;

    template<class S>
    bool feasible(S const& solution, double tolerance=1e-3) const;

private:
    ConstraintSense _sense;
    E _expression;
};

}

#include "RawConstraint.tplt"

#endif //ROBUSTOPTIMIZATION_RAWCONSTRAINT_H
