//
// Created by simon on 14/06/2021.
//

#ifndef ROBUSTOPTIMIZATION_OBJECTIVEBASE_H
#define ROBUSTOPTIMIZATION_OBJECTIVEBASE_H

#include <string>

#include "predeclarations_and_concepts.h"

namespace robust_model {

template<class E>
class ObjectiveBase {
public:
    using Sense = ObjectiveSense;

public:
    ObjectiveBase() = default;

    ObjectiveBase(Sense sense, E const& expression);

    std::string to_string() const;

    Sense sense() const;

    E const& expression() const;

    double value() const;

    template<class S>
    double value(S const& s) const;

    void set_solution(double solution);
    void invalidate_solution();

private:
    Sense _sense;
    E _expression;
    std::optional<double> _value;
};

}

#include "Objective.tplt"

#endif //ROBUSTOPTIMIZATION_OBJECTIVEBASE_H
