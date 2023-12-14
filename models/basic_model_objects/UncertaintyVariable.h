//
// Created by simon on 14/06/2021.
//

#ifndef ROBUSTOPTIMIZATION_UNCERTAINTYVARIABLE_H
#define ROBUSTOPTIMIZATION_UNCERTAINTYVARIABLE_H

#include <optional>

#include "types_and_constants.h"
#include "VariableBase.h"

namespace robust_model {

class UncertaintyVariable : public VariableBase<UncertaintyVariable> {
public:
    UncertaintyVariable(Index id, std::string const& name,
                        std::optional<period_id> p, double lb,
                        double ub);

    bool has_period() const;

    period_id period() const;

    bool is_0_1_var() const;

private:
    std::optional<period_id> const _period;
};

}

#endif //ROBUSTOPTIMIZATION_UNCERTAINTYVARIABLE_H
