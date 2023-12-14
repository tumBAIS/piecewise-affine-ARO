//
// Created by simon on 11/08/2021.
//

#include "ParallelInstanceEvaluator.h"

namespace testing {
ParallelInstanceEvaluator::ParallelInstanceEvaluator(std::ostream& output_stream, InstanceGeneratorBase& generator) :
        _output_stream(output_stream), _instance_generator(generator) {}

void ParallelInstanceEvaluator::set_num_simulations(size_t num_simulations) {
    _num_simulations = num_simulations;
}

void ParallelInstanceEvaluator::add_test(std::string const& name, std::function<std::tuple<double, double>(
        robust_model::ROModel const&)> const& test) {
    _tests.emplace_back(test, name);
}

void ParallelInstanceEvaluator::add_simulation_test(std::string const& name,
                                                    std::function<std::tuple<double, double, double>(
                                                            robust_model::ROModel const&,
                                                            std::vector<std::vector<double>> const&)> const& test) {
    _simulation_tests.emplace_back(test, name);
}

void ParallelInstanceEvaluator::thread_function() {
    while (true) {
        auto const instance = _instance_generator.next_instance();
        if (not instance.first) {
            return;
        }
        for (auto const& [test, test_description]: _tests) {
            auto const result = test(*instance.first);
            _output_lock.lock();
            _output_stream << instance.second << test_description << ";" << std::get<0>(result) << ";"
                           << std::get<1>(result) << ";" << std::nan("0") << std::endl;
            _output_lock.unlock();
        }
        if (_num_simulations.has_value()) {
            auto const simulated_scenarios = _simulated_realization_generator(instance.first->uncertainty_set(),
                                                                              _num_simulations.value());
            for (auto const& [test, test_description]: _simulation_tests) {
                auto const result = test(*instance.first, simulated_scenarios);
                _output_lock.lock();
                _output_stream << instance.second << test_description << ";" << std::get<0>(result) << ";"
                               << std::get<1>(result) << ";" << std::get<2>(result) << std::endl;
                _output_lock.unlock();
            }
        } else {
            helpers::exception_check(_simulation_tests.empty(),
                                     "Existing simulation tests with no number of simulated runs");
        }
    }
}

void ParallelInstanceEvaluator::run_tests(size_t num_threads) {
    auto t = std::chrono::system_clock::now();
    _output_stream << _instance_generator.descriptions() << "test_name;runtime;objective;avg_samples" << std::endl;
    std::vector<std::thread> threads;
    for (size_t i = 0; i < num_threads - 1; ++i) {
        threads.emplace_back(&ParallelInstanceEvaluator::thread_function, this);
    }
    thread_function();
    for (auto& thread: threads) {
        thread.join();
    }
    helpers::global_logger << "Test time: " + std::to_string(
            double(std::chrono::duration_cast<std::chrono::milliseconds>(
                    std::chrono::system_clock::now() - t).count()) / 1000) + "s";
}

void ParallelInstanceEvaluator::set_simulated_realization_generator(
        std::function<std::vector<std::vector<double>>(robust_model::UncertaintySet const&, size_t)> const& generator) {
    _simulated_realization_generator = generator;
}


}
