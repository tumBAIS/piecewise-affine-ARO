//
// Created by simon on 11/08/2021.
//

#include "InstanceGeneratorBase.h"

namespace testing {
void InstanceGeneratorBase::add_set_type(robust_model::UncertaintySet::SpecialSetType t) {
    _set_types.emplace_back(t);
}

void InstanceGeneratorBase::add_budget(std::function<double(size_t)> const& generator, std::string const& description) {
    _budgets.emplace_back(generator, description);
}

void InstanceGeneratorBase::set_number_of_iterations(size_t number_of_iterations) {
    _number_of_iterations = number_of_iterations;
}

std::pair<std::unique_ptr<robust_model::ROModel>, std::string> InstanceGeneratorBase::next_instance() {
    _generation_lock.lock();
    if (_active) {
        std::pair<std::unique_ptr<robust_model::ROModel>, std::string>
                instance{generate_instance(), instance_description()};
        auto const b = budget(instance.first->num_uvars());
        if(b.has_value()){
            instance.first->add_special_uncertainty_constraint(set_type(), b.value());
        }
        else{
            instance.first->add_special_uncertainty_constraint(set_type());
        }
        _active = increment() and _active;
        _generation_lock.unlock();
        return instance;
    } else {
        _generation_lock.unlock();
        return {};
    }
}

std::string InstanceGeneratorBase::descriptions() const {
    return "uncertainty_type;budget;" + descriptions_test_specific() + "iteration;";
}

std::string InstanceGeneratorBase::instance_description() {
    std::string s;
    s += robust_model::UncertaintySet::to_string(_set_types.at(_set_types_id)) + ";";
    s += budget_description() + ";";
    s += instance_description_test_specific();
    s += std::to_string(_iteration) + ";";
    return s;
}

bool InstanceGeneratorBase::increment() {

    if (++_iteration < _number_of_iterations) {
        return true;
    }
    _iteration = 0;
    helpers::global_logger << "Finished instances in block " + instance_description();
    if (increment_test_specific()) {
        return true;
    }
    if (_set_types.at(_set_types_id) == robust_model::UncertaintySet::SpecialSetType::BUDGET) {
        if (++_budgets_id < _budgets.size()) {
            return true;
        }
        _budgets_id = 0;
    }
    if (++_set_types_id < _set_types.size()) {
        return true;
    }
    _set_types_id = 0;
    return false;
}

robust_model::UncertaintySet::SpecialSetType InstanceGeneratorBase::set_type() const {
    return _set_types.at(_set_types_id);
}

std::string InstanceGeneratorBase::budget_description() const {
    return (set_type() == robust_model::UncertaintySet::SpecialSetType::BUDGET) ? _budgets.at(_budgets_id).second : "";
}

std::optional<double> InstanceGeneratorBase::budget(size_t m) const {
    return (set_type() == robust_model::UncertaintySet::SpecialSetType::BUDGET) ?
           std::optional<double>{_budgets.at(_budgets_id).first(m)} : std::optional<double>{};
}

}