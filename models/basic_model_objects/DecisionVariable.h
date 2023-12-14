//
// Created by simon on 14/06/2021.
//

#ifndef ROBUSTOPTIMIZATION_DECISIONVARIABLE_H
#define ROBUSTOPTIMIZATION_DECISIONVARIABLE_H

#include <optional>
#include <map>

#include "../../helpers/helpers.h"
#include "types_and_constants.h"
#include "VariableBase.h"
#include "UncertaintyVariable.h"

namespace robust_model {

struct DominatingSolution {
    explicit DominatingSolution(std::vector<double> solutions);

    std::vector<double> const& values_per_vertex() const;

    std::string to_string() const;

private:
    std::vector<double> _values_per_vertex;
};

struct AffineSolution {
    explicit AffineSolution(double constant,
                            std::map<UncertaintyVariable::Index, double> dependent_scales);

    double constant() const;

    std::map<UncertaintyVariable::Index, double> const& dependent_scales() const;

    std::string to_string() const;

private:
    double _constant;
    std::map<UncertaintyVariable::Index, double> _dependent_scales;
};

class DecisionVariable : public VariableBase<DecisionVariable> {
public:
    DecisionVariable(Index id, std::string const& name,
                     std::optional<period_id> p,
                     double lb, double ub);

    void add_dependency(UncertaintyVariable::Index dependency);

    void add_dependencies(std::vector<UncertaintyVariable::Reference> const& dependencies);

    std::vector<UncertaintyVariable::Reference> const& dependencies() const;

    bool has_period() const;

    period_id period() const;

    AffineSolution const& affine_solution() const;

    void set_affine_solution(AffineSolution const& affine_solution);

    DominatingSolution const& dominating_solution() const;

    void set_dominating_solution(DominatingSolution const& dominating_solution);

private:
    std::vector<UncertaintyVariable::Reference> _dependencies;
    std::optional<period_id> const _period;
    std::optional<AffineSolution> _affine_solution;
    std::optional<DominatingSolution> _dominating_solution;
};

}

#endif //ROBUSTOPTIMIZATION_DECISIONVARIABLE_H
