//
// Created by simon on 11/08/2021.
//

#ifndef ROBUSTOPTIMIZATION_INSTANCEGENERATORBASE_H
#define ROBUSTOPTIMIZATION_INSTANCEGENERATORBASE_H

#include "../../models/ROModel.h"
#include <mutex>
#include <memory>
#include <functional>

namespace testing {

class InstanceGeneratorBase {
public:
    void add_set_type(robust_model::UncertaintySet::SpecialSetType t);

    void add_budget(std::function<double(size_t)> const& generator, std::string const& description);

    void set_number_of_iterations(size_t number_of_iterations);

    std::pair<std::unique_ptr<robust_model::ROModel>, std::string> next_instance();

    std::string descriptions() const;

private:
    virtual std::string descriptions_test_specific() const = 0;

    virtual std::unique_ptr<robust_model::ROModel> generate_instance() = 0;

    std::string instance_description();
    virtual std::string instance_description_test_specific() = 0;

    bool increment();
    virtual bool increment_test_specific() = 0;

protected:
    robust_model::UncertaintySet::SpecialSetType set_type() const;

    std::string budget_description() const;
    std::optional<double> budget(size_t) const;

private:
    bool _active = true;
    std::mutex _generation_lock;
    std::vector<robust_model::UncertaintySet::SpecialSetType> _set_types;
    std::vector<std::pair<std::function<double(size_t)>, std::string>> _budgets;
    size_t _number_of_iterations;

    size_t _set_types_id = 0,
            _budgets_id = 0;
    size_t _iteration = 0;
};

}

#endif //ROBUSTOPTIMIZATION_INSTANCEGENERATORBASE_H
