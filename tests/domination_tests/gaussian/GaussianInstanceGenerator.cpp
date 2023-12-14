//
// Created by simon on 11/08/2021.
//

#include "GaussianInstanceGenerator.h"

#include <chrono>
#include <random>

namespace testing {

std::unique_ptr<robust_model::ROModel> GaussianInstanceGenerator::generate_instance() {
    auto m = _test_sizes.at(_test_sizes_id);
    auto vars_per_period = _num_vars_per_period.at(_num_vars_per_period_id).first(m);
    auto objective_scale = _objective_uncertainty_scales.at(_objective_uncertainty_scales_id);

    std::ranlux48 generator(std::chrono::system_clock::now().time_since_epoch().count());
    std::normal_distribution<double> distribution(0, 1);

    auto model_uptr = std::make_unique<robust_model::ROModel>("MyModel");
    auto& model = *model_uptr;
    auto uvs = model.add_uncertainty_variables(0, "UVar", 0, 0, 1);
    auto dvs = model.add_decision_variables(0, "Var", 0, 0);
    for (size_t i = 0; i < m; ++i) {
        uvs.emplace_back(model.add_uncertainty_variable("UVar" + std::to_string(i), i / vars_per_period, 0., 1.));
        dvs.emplace_back(model.add_decision_variable("DVar" + std::to_string(i), i / vars_per_period, 0.));
    }
    robust_model::RoAffineExpression objective;
    for (int i = 0; i < m; ++i) {
        objective += (1 + objective_scale * std::abs(distribution(generator))) * dvs[i];
        robust_model::LinearExpression<robust_model::DecisionVariable::Reference> lhs;
        robust_model::LinearExpression<robust_model::UncertaintyVariable::Reference> rhs;
        for (int j = 0; j < m; ++j) {
            lhs += (int(i == j) + std::abs(distribution(generator)) / std::sqrt(m)) * dvs[j];
//            rhs += (int(i == j) + std::abs(distribution(generator)) / std::sqrt(m)) * uvs[j];
        }
        model.add_constraint(lhs >= uvs[i], "C" + std::to_string(i));
    }
    model.set_objective(objective, robust_model::ObjectiveSense::MIN);

    return model_uptr;
}

std::string GaussianInstanceGenerator::instance_description_test_specific() {
    std::string s;
    s += _num_vars_per_period.at(_num_vars_per_period_id).second + ";";
    s += std::to_string(_objective_uncertainty_scales.at(_objective_uncertainty_scales_id)) + ";";
    s += std::to_string(_test_sizes.at(_test_sizes_id)) + ";";
    return s;
}

bool GaussianInstanceGenerator::increment_test_specific() {
    {
        if (++_test_sizes_id < _test_sizes.size()) {
            return true;
        }
        _test_sizes_id = 0;
        if (++_objective_uncertainty_scales_id < _objective_uncertainty_scales.size()) {
            return true;
        }
        _objective_uncertainty_scales_id = 0;
        if (++_num_vars_per_period_id < _num_vars_per_period.size()) {
            return true;
        }
        _num_vars_per_period_id = 0;
        return false;
    }
}

void GaussianInstanceGenerator::add_test_size(size_t m) {
    _test_sizes.emplace_back(m);
}

void GaussianInstanceGenerator::add_num_vars_per_period_generator(std::function<size_t(size_t)> const& generator,
                                                                  std::string const& description) {
    _num_vars_per_period.emplace_back(generator, description);
}

void GaussianInstanceGenerator::add_objective_uncertainty_scale(double s) {
    _objective_uncertainty_scales.emplace_back(s);
}


std::string GaussianInstanceGenerator::descriptions_test_specific() const {
    return "variables_per_period;objective_uncertainty_scale;test_size;";
}

}