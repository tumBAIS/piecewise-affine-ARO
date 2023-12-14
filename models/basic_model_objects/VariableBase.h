//
// Created by simon on 14/06/2021.
//

#ifndef ROBUSTOPTIMIZATION_VARIABLEBASE_H
#define ROBUSTOPTIMIZATION_VARIABLEBASE_H

#include "../../helpers/helpers.h"
#include "types_and_constants.h"
#include "VariableReference.h"

namespace robust_model {

template<class V>
class VariableBase : public helpers::IndexedObject<V> {
public:
    using Index = typename helpers::IndexedObject<V>::Index;
    using Reference = VariableReference<V>;

public:
    VariableBase(Index id, std::string name, double lb, double ub);

    Reference reference() const;

    std::string const& name() const;

    double lb() const;

    double ub() const;

    bool bounded() const;

private:
    std::string const _name;
    double const _lb;
    double const _ub;
};

}

#include "VariableBase.tplt"

#endif //ROBUSTOPTIMIZATION_VARIABLEBASE_H
