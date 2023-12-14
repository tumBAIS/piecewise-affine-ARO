//
// Created by simon on 06/07/2021.
//

#include "gaussian/GaussianInstanceGenerator.h"
#include "../test_helpers/ParallelInstanceEvaluator.h"
#include "../../solvers/aro_policy_solvers/DominatingPolicySolver.h"
#include "../../solvers/aro_policy_solvers/AffineAdjustablePolicySolver.h"
#include "../../solvers/aro_policy_solvers/LiftingPolicySolver.h"
#include "../../solvers/aro_policy_solvers/BGPiecewiseLinearPolicySolver.h"
#include "../../helpers/Logger.h"

#include <iostream>
#include <fstream>
#include <cmath>

int
main(int argc,
     char *argv[]) {
//    try{
    std::string run_name = "gaussian_test_" + helpers::time_stamp();
    helpers::global_logger.set_logfile("logs/" + run_name + ".log");
    helpers::global_logger << "Logging " + run_name;
    std::ofstream output_stream("results/" + run_name + ".csv");
    auto instance_generator = testing::GaussianInstanceGenerator();
    auto tester = testing::ParallelInstanceEvaluator(output_stream, instance_generator);

    double const time_limit = 14400;

    tester.add_test("PAPBT",
                    [time_limit](robust_model::ROModel const& model) {
                        auto dm = robust_model::DominatingPolicySolver(
                                model, robust_model::DominatingPolicySolver::Mode::BENTAL);
                        dm.build();
                        dm.set_runtime_limit(time_limit);
                        dm.solve();
                        return std::make_tuple(dm.runtime(), dm.objective_value());
                    });
    tester.add_test("PAP",
                    [time_limit](robust_model::ROModel const& model) {
                        auto dm = robust_model::DominatingPolicySolver(
                                model, robust_model::DominatingPolicySolver::Mode::THOMAE);
                        dm.build();
                        dm.set_runtime_limit(time_limit);
                        dm.solve();
                        return std::make_tuple(dm.runtime(), dm.objective_value());
                    });
    tester.add_test("SPAP",
                    [time_limit](robust_model::ROModel const& model) {
                        auto dm = robust_model::DominatingPolicySolver(
                                model, robust_model::DominatingPolicySolver::Mode::THOMAE_SHIFT);
                        dm.build();
                        dm.set_runtime_limit(time_limit);
                        dm.solve();
                        return std::make_tuple(dm.runtime(), dm.objective_value());
                    });
    tester.add_test("AFF",
                    [time_limit](robust_model::ROModel const& model) {
                        auto am = robust_model::AffineAdjustablePolicySolver(model);
                        am.build();
                        am.set_runtime_limit(time_limit);
                        am.solve();
                        return std::make_tuple(am.runtime(), am.objective_value());
                    });
    tester.add_test("TLIFT",
                    [time_limit](robust_model::ROModel const& model) {
                        auto lm = robust_model::LiftingPolicySolver(model);
                        lm.add_domination_motivated_break_points();
                        lm.build();
                        lm.set_runtime_limit(time_limit);
                        lm.solve();
                        return std::make_tuple(lm.runtime(), lm.objective_value());
                    });
    tester.add_test("BG",
                    [time_limit](robust_model::ROModel const& model) {
                        auto bgm = robust_model::BGPiecewiseLinearPolicySolver(model, 2);
                        bgm.build();
                        bgm.set_runtime_limit(time_limit);
                        bgm.solve();
                        return std::make_tuple(bgm.runtime(), bgm.objective_value());
                    });
    tester.add_test("BOX",
                    [time_limit](robust_model::ROModel const& model) {
                        auto dm = robust_model::DominatingPolicySolver(
                                model, robust_model::DominatingPolicySolver::Mode::BOX);
                        dm.build();
                        dm.set_runtime_limit(time_limit);
                        dm.solve();
                        return std::make_tuple(dm.runtime(), dm.objective_value());
                    });

    for (size_t i = 2; i <= 10; ++i) {
        instance_generator.add_test_size(i * i);
    }

    instance_generator.add_set_type(robust_model::UncertaintySet::SpecialSetType::BALL);
    instance_generator.add_set_type(robust_model::UncertaintySet::SpecialSetType::BUDGET);

    instance_generator.add_objective_uncertainty_scale(0);
    instance_generator.add_objective_uncertainty_scale(.1);
    instance_generator.add_objective_uncertainty_scale(.5);
    instance_generator.add_objective_uncertainty_scale(1);
    instance_generator.add_objective_uncertainty_scale(5);

    instance_generator.add_budget([](size_t m) { return std::sqrt(m); }, "sqrt(m)");

    instance_generator.add_num_vars_per_period_generator([](size_t m) { return size_t(std::ceil(std::sqrt(m))); },
                                                         "sqrt(m)");

    instance_generator.set_number_of_iterations(30);

    tester.run_tests(6);
//    }
//    catch(std::exception &e){
//        helpers::exception_check(false, e.what());
//    }
}