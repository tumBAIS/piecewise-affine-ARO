//
// Created by simon on 25/10/2021.
//

#include "DemandCoveringInstanceGenerator.h"
#include "DemandCoveringProblem.h"

#include <random>
#include <chrono>

namespace testing{

void DemandCoveringInstanceGenerator::add_num_locations(size_t m) {
    _numbers_of_locations.emplace_back(m);
}

void DemandCoveringInstanceGenerator::add_num_periods(size_t m) {
    _numbers_of_periods.emplace_back(m);
}

void DemandCoveringInstanceGenerator::add_lost_demand_cost(double c) {
    _lost_demand_costs.emplace_back(c);
}

std::string DemandCoveringInstanceGenerator::descriptions_test_specific() const {
    return "num_locations;num_periods;lost_demand_cost;";
}

std::unique_ptr<robust_model::ROModel> DemandCoveringInstanceGenerator::generate_instance() {
    auto const num_periods = _numbers_of_periods.at(_numbers_of_periods_id);
    auto const lost_demands_cost = _lost_demand_costs.at(_lost_demand_costs_id);
    auto const num_locations = _numbers_of_locations.at(_numbers_of_locations_id);

    demand_covering_problem::DemandCoveringProblem problem(num_periods,
                                                           8,
                                                           1,
                                                           0,
                                                           .02,
                                                           .1,
                                                           .2,
                                                           lost_demands_cost);

    std::ranlux48 generator(std::chrono::system_clock::now().time_since_epoch().count());
    std::uniform_int_distribution<> uniform_distribution(0, 2*int(sqrt(num_locations))+1);

    for(size_t i = 0; i < num_locations; ++i){
        problem.add_location(uniform_distribution(generator), uniform_distribution(generator));
    }

    problem.set_demand_pattern([](auto period, auto location){
        std::ranlux48 generator(std::chrono::system_clock::now().time_since_epoch().count());
        std::normal_distribution<double> normal_distribution(0, 1);

        return 10 + 2* normal_distribution(generator);
    });

    demand_covering_problem::DemandCoveringROModelBuilder math_model_builder(problem);

    math_model_builder.build();
    auto model = math_model_builder.model_unique_ptr();
    return model;
}

std::string DemandCoveringInstanceGenerator::instance_description_test_specific() {
    std::string s;
    s += std::to_string(_numbers_of_locations.at(_numbers_of_locations_id)) + ";";
    s += std::to_string(_numbers_of_periods.at(_numbers_of_periods_id)) + ";";
    s += std::to_string(_lost_demand_costs.at(_lost_demand_costs_id)) + ";";
    return s;}

bool DemandCoveringInstanceGenerator::increment_test_specific() {
    if (++_numbers_of_periods_id < _numbers_of_periods.size()) {
        return true;
    }
    _numbers_of_periods_id = 0;
    if (++_numbers_of_locations_id < _numbers_of_locations.size()) {
        return true;
    }
    _numbers_of_locations_id = 0;
    if (++_lost_demand_costs_id < _lost_demand_costs.size()) {
        return true;
    }
    _lost_demand_costs_id = 0;
    return false;}

}
