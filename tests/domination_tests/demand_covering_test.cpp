//
// Created by simon on 28/05/2021.
//

#include "iostream"
#include "random"

#include "../../solvers/aro_policy_solvers/DominatingPolicySolver.h"
#include "../../solvers/aro_policy_solvers/AffineAdjustablePolicySolver.h"
#include "../../solvers/aro_policy_solvers/LiftingPolicySolver.h"
#include "../../solvers/aro_policy_solvers/BGPiecewiseLinearPolicySolver.h"
#include "../test_helpers/ParallelInstanceEvaluator.h"
#include "demand_covering/DemandCoveringInstanceGenerator.h"

int
main(int argc,
     char *argv[]) {
//    try {
    std::string run_name = "demand_covering_test" + helpers::time_stamp();
    helpers::global_logger.set_logfile("logs/" + run_name + ".log");
    helpers::global_logger << "Logging " + run_name;
    std::ofstream output_stream("results/" + run_name + ".csv");
    auto instance_generator = testing::DemandCoveringInstanceGenerator();
    for (size_t i = 1; i <= 5; ++i) {
        instance_generator.add_num_locations(i * 2);
    }

    for (size_t i = 0; i <= 3; ++i) {
        instance_generator.add_num_periods(2 * i + 1);
    }

    instance_generator.add_lost_demand_cost(.1);
    instance_generator.add_lost_demand_cost(.25);
    instance_generator.add_lost_demand_cost(.5);

    auto tester = testing::ParallelInstanceEvaluator(output_stream, instance_generator);

    double const time_limit = 7200;

    tester.add_simulation_test("PAPBT",
                               [time_limit](robust_model::ROModel const& model,
                                  std::vector<std::vector<double>> const& realizations) {
                                   auto dm = robust_model::DominatingPolicySolver(
                                           model, robust_model::DominatingPolicySolver::Mode::BENTAL);
                                   dm.set_reoptimize(true);
                                   dm.build();
                                   dm.set_runtime_limit(time_limit);
                                   dm.solve();
                                   double ov_avg = 0;
                                   for (auto const& realization: realizations) {
                                       ov_avg += model.objective_value_for_solution(
                                               dm.specific_solution(realization));
                                   }
                                   return std::make_tuple(dm.runtime(), dm.objective_value(),
                                                          ov_avg / double(realizations.size()));
                               });
    tester.add_simulation_test("PAP",
                               [time_limit](robust_model::ROModel const& model,
                                  std::vector<std::vector<double>> const& realizations) {
                                   auto dm = robust_model::DominatingPolicySolver(
                                           model, robust_model::DominatingPolicySolver::Mode::THOMAE);
                                   dm.set_reoptimize(true);
                                   dm.build();
                                   dm.set_runtime_limit(time_limit);
                                   dm.solve();
                                   double ov_avg = 0;
                                   for (auto const& realization: realizations) {
                                       ov_avg += model.objective_value_for_solution(
                                               dm.specific_solution(realization));
                                   }
                                   return std::make_tuple(dm.runtime(), dm.objective_value(),
                                                          ov_avg / double(realizations.size()));
                               });
    tester.add_simulation_test("SPAP",
                               [time_limit](robust_model::ROModel const& model,
                                  std::vector<std::vector<double>> const& realizations) {
                                   auto dm = robust_model::DominatingPolicySolver(
                                           model, robust_model::DominatingPolicySolver::Mode::THOMAE_SHIFT);
                                   dm.set_reoptimize(true);
                                   dm.build();
                                   dm.set_runtime_limit(time_limit);
                                   dm.solve();
                                   double ov_avg = 0;
                                   for (auto const& realization: realizations) {
                                       ov_avg += model.objective_value_for_solution(
                                               dm.specific_solution(realization));
                                   }
                                   return std::make_tuple(dm.runtime(), dm.objective_value(),
                                                          ov_avg / double(realizations.size()));
                               });
    tester.add_simulation_test("AFF",
                               [time_limit](robust_model::ROModel const& model,
                                  std::vector<std::vector<double>> const& realizations) {
                                   auto am = robust_model::AffineAdjustablePolicySolver(model);
                                   am.build();
                                   am.set_runtime_limit(time_limit);
                                   am.add_average_reoptimization();
                                   am.solve();
                                   double ov_avg = 0;
                                   for (auto const& realization: realizations) {
                                       ov_avg += model.objective_value_for_solution(
                                               am.specific_solution(realization));
                                   }
                                   return std::make_tuple(am.runtime(), am.objective_value(),
                                                          ov_avg / double(realizations.size()));
                               });
    tester.add_simulation_test("TLIFT",
                               [time_limit](robust_model::ROModel const& model,
                                  std::vector<std::vector<double>> const& realizations) {
                                   auto lm = robust_model::LiftingPolicySolver(model);
                                   lm.add_domination_motivated_break_points();
                                   lm.build();
                                   lm.set_runtime_limit(time_limit);
                                   lm.add_average_lifted_vertex_reoptimization();
                                   lm.solve();
                                   double ov_avg = 0;
                                   for (auto const& realization: realizations) {
                                       ov_avg += model.objective_value_for_solution(
                                               lm.specific_solution(realization));
                                   }
                                   return std::make_tuple(lm.runtime(), lm.objective_value(),
                                                          ov_avg / double(realizations.size()));
                               });
    tester.add_simulation_test("BG",
                               [time_limit](robust_model::ROModel const& model,
                                  std::vector<std::vector<double>> const& realizations) {
                                   auto bgm = robust_model::BGPiecewiseLinearPolicySolver(model, 2);
                                   bgm.build();
                                   bgm.set_runtime_limit(time_limit);
                                   bgm.solve();
                                   double ov_avg = 0;
                                   for (auto const& realization: realizations) {
                                       ov_avg += model.objective_value_for_solution(
                                               bgm.specific_solution(realization));
                                   }
                                   return std::make_tuple(bgm.runtime(), bgm.objective_value(),
                                                          ov_avg / double(realizations.size()));
                               });
    tester.add_simulation_test("BOX",
                               [time_limit](robust_model::ROModel const& model,
                                  std::vector<std::vector<double>> const& realizations) {
                                   auto dm = robust_model::DominatingPolicySolver(
                                           model, robust_model::DominatingPolicySolver::Mode::BOX);
                                   dm.build();
                                   dm.set_runtime_limit(time_limit);
                                   dm.solve();
                                   double ov_avg = 0;
                                   for (auto const& realization: realizations) {
                                       ov_avg += model.objective_value_for_solution(
                                               dm.specific_solution(realization));
                                   }
                                   return std::make_tuple(dm.runtime(), dm.objective_value(),
                                                          ov_avg / double(realizations.size()));
                               });

    instance_generator.add_set_type(robust_model::UncertaintySet::SpecialSetType::BALL);
    instance_generator.add_set_type(robust_model::UncertaintySet::SpecialSetType::BUDGET);

    instance_generator.add_budget([](size_t m) { return std::sqrt(m); }, "sqrt(m)");

    instance_generator.set_number_of_iterations(25);
    tester.set_num_simulations(500);

    tester.run_tests(6);
//    }
//    catch (std::exception& e) {
//        helpers::exception_check(false, e.what());
//    }
}