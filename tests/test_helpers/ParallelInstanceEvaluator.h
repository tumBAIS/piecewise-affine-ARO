//
// Created by simon on 11/08/2021.
//

#ifndef ROBUSTOPTIMIZATION_PARALLELINSTANCEEVALUATOR_H
#define ROBUSTOPTIMIZATION_PARALLELINSTANCEEVALUATOR_H

#include "InstanceGeneratorBase.h"
#include <functional>
#include <thread>
#include <cmath>

namespace testing {

class ParallelInstanceEvaluator {
public:
    explicit ParallelInstanceEvaluator(std::ostream& output_stream, InstanceGeneratorBase& generator);

    void set_num_simulations(size_t num_simulations);

    void add_test(std::string const& name,
                  std::function<std::tuple<double, double>(robust_model::ROModel const&)> const& test);

    void add_simulation_test(
            std::string const& name,
            std::function<std::tuple<double, double, double>(robust_model::ROModel const&,
                                                             std::vector<std::vector<double>> const&)> const& test);

    void thread_function();

    void run_tests(size_t num_threads = 1);

    void set_simulated_realization_generator(
            std::function<std::vector<std::vector<double>>(robust_model::UncertaintySet const&, size_t)> const&
            generator);

private:
    InstanceGeneratorBase& _instance_generator;
    std::vector<std::pair<
            std::function<std::tuple<double, double>(robust_model::ROModel const&)>,
            std::string>> _tests;
    std::vector<std::pair<
            std::function<std::tuple<double, double, double>(robust_model::ROModel const&,
                                                             std::vector<std::vector<double>> const&)>,
            std::string>> _simulation_tests;
    std::mutex _output_lock;
    std::ostream& _output_stream;
    std::optional<size_t> _num_simulations = std::optional<size_t>{};
    std::function<std::vector<std::vector<double>>(robust_model::UncertaintySet const&, size_t)>
            _simulated_realization_generator = [](robust_model::UncertaintySet const& uncertainty_set,
                                                  size_t num_realizations) {
        return uncertainty_set.generate_uncertainty(num_realizations);
    };
};

}

#endif //ROBUSTOPTIMIZATION_PARALLELINSTANCEEVALUATOR_H
