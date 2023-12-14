//
// Created by simon on 14/06/2021.
//

#ifndef ROBUSTOPTIMIZATION_SOCVARIABLE_H
#define ROBUSTOPTIMIZATION_SOCVARIABLE_H

#include <optional>

#include "VariableBase.h"
#include "types_and_constants.h"

namespace robust_model {

class SOCVariable;

class DualInformation {
public:
    enum class DualType{
        CONSTRAINT,
        CONSTRAINT_NORMED_VECTOR,
        VARIABLE_LB,
        VARIABLE_UB
    };
public:
    DualInformation(DualType type, size_t constraint_number);
    DualInformation(DualType type, size_t constraint_number, size_t normed_vector_entry);
    DualInformation(DualType type, helpers::SmartIndex<SOCVariable> variable);

    bool is_constraint_dual() const;
    bool is_normed_vector_dual() const;
    bool is_lower_bound_dual() const;
    bool is_upper_bound_dual() const;

    size_t constraint_number() const;
    size_t normed_vector_entry() const;
    helpers::SmartIndex<SOCVariable> bounded_dual_variable() const;

private:
    DualType const _type;
    std::optional<size_t> const _constraint_number;
    std::optional<size_t> const _normed_vector_entry;
    std::optional<helpers::SmartIndex<SOCVariable>> const _variable;
};

class SOCVariable : public VariableBase<SOCVariable> {
public:
    SOCVariable(Index id, std::string const& name, double lb, double ub, VariableType type=VariableType::Continuous);

    SOCVariable(Index id, std::string const& name, double lb, double ub, DualInformation const& dual_info);

    SOCVariable(Index id, std::string const& name, double lb, double ub, std::optional<DualInformation> const& dual_info,
                VariableType type);

    void set_solution(double sol);
    void invalidate_solution();

    bool has_solution() const;

    double solution() const;

    VariableType const type() const;

    bool has_dual_info() const;
    DualInformation const& dual_information() const;

private:
    std::optional<double> _solution;
    VariableType const _type;
    std::optional<DualInformation> const _dual_info;
};

}


#endif //ROBUSTOPTIMIZATION_SOCVARIABLE_H
