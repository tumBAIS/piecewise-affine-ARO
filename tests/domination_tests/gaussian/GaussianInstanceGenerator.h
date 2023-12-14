//
// Created by simon on 11/08/2021.
//

#ifndef ROBUSTOPTIMIZATION_GAUSSIANINSTANCEGENERATOR_H
#define ROBUSTOPTIMIZATION_GAUSSIANINSTANCEGENERATOR_H

#include "../../test_helpers/InstanceGeneratorBase.h"
#include <functional>

namespace testing {

class GaussianInstanceGenerator : public InstanceGeneratorBase {
public:
    void add_test_size(size_t m);
    void add_num_vars_per_period_generator(std::function<size_t(size_t)> const& generator,
                                           std::string const& description);

    void add_objective_uncertainty_scale(double s);

private:
    std::string descriptions_test_specific() const override;

    std::unique_ptr<robust_model::ROModel> generate_instance() override;

    std::string instance_description_test_specific() override;

    bool increment_test_specific() override;

private:
    std::vector<std::pair<std::function<size_t(size_t)>, std::string>> _num_vars_per_period;
    std::vector<double> _objective_uncertainty_scales;
    std::vector<size_t> _test_sizes;

    size_t _num_vars_per_period_id = 0,
            _objective_uncertainty_scales_id = 0,
            _test_sizes_id = 0;
};

}


#endif //ROBUSTOPTIMIZATION_GAUSSIANINSTANCEGENERATOR_H
