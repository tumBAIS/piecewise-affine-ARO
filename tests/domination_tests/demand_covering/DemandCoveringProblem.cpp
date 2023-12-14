//
// Created by simon on 20/10/2021.
//

#include "DemandCoveringProblem.h"
#include <algorithm>
#include <cmath>

namespace demand_covering_problem {


Location::Location(helpers::SmartIndex<Location> const& id, double x_coord, double y_coord) : IndexedObject(id),
                                                                                              _x_coord(x_coord),
                                                                                              _y_coord(y_coord) {}

double Location::distance_to(Location::Index other) const {
    return std::sqrt(std::pow(x_coord() - other->x_coord(), 2) +
                std::pow(y_coord() - other->y_coord(), 2));
}

double Location::x_coord() const {
    return _x_coord;
}

double Location::y_coord() const {
    return _y_coord;
}

DemandCoveringProblem::DemandCoveringProblem(size_t const num_planning_periods,
                                             size_t const num_execution_per_planning_periods, double const capacity_cost,
                                             double const capacity_usage_cost, double const distance_loss_factor,
                                             double const next_period_loss, double const next_day_loss,
                                             double const lost_demand_cost) : _num_planning_periods(num_planning_periods),
                                                                        _num_execution_per_planning_periods(
                                                                                num_execution_per_planning_periods),
                                                                        _capacity_cost(capacity_cost),
                                                                        _capacity_usage_cost(capacity_usage_cost),
                                                                        _distance_loss_factor(distance_loss_factor),
                                                                        _next_period_loss(next_period_loss),
                                                                        _next_day_loss(next_day_loss),
                                                                        _lost_demand_cost(lost_demand_cost),
                                                                        _base_demand_pattern(
                                                                                [](robust_model::period_id , Location::Index) { return 1; }) {}

helpers::SmartIndex<Location> DemandCoveringProblem::add_location(double x, double y) {
    return IndexedObjectOwner<Location>::base_add_object(x, y);
}

void
DemandCoveringProblem::set_demand_pattern(std::function<double(robust_model::period_id, Location::Index)> demand_pattern) {
    _base_demand_pattern = std::move(demand_pattern);
}

double DemandCoveringProblem::base_demand(robust_model::period_id period, Location::Index location) const {
    return _base_demand_pattern(period, location);
}

size_t DemandCoveringProblem::num_locations() const {
    return IndexedObjectOwner<Location>::num_objects();
}

std::vector<Location> const& DemandCoveringProblem::locations() const {
    return IndexedObjectOwner<Location>::objects();
}

std::vector<Location::Index> const& DemandCoveringProblem::location_ids() const {
    return IndexedObjectOwner<Location>::ids();
}

size_t DemandCoveringProblem::num_planning_periods() const {
    return _num_planning_periods;
}

size_t DemandCoveringProblem::num_execution_periods_per_planning_period() const {
    return _num_execution_per_planning_periods;
}

size_t DemandCoveringProblem::num_periods_per_planning_period() const {
    return num_execution_periods_per_planning_period() + 1;
}

double DemandCoveringProblem::capacity_cost() const {
    return _capacity_cost;
}

double DemandCoveringProblem::capacity_usage_cost() const {
    return _capacity_usage_cost;
}

double DemandCoveringProblem::distance_loss_factor(double const dist) const {
    return std::min(_distance_loss_factor * dist, 1.);
}

double DemandCoveringProblem::next_period_loss_factor(robust_model::period_id const period) const {
    if (period + 1 == planning_period_period_id(period_to_execution_period(period + 1))) {
        return _next_day_loss;
    } else {
        return _next_period_loss;
    }
}

robust_model::period_id DemandCoveringProblem::planning_period_period_id(size_t day) const {
    return robust_model::period_id(day * (num_periods_per_planning_period()));
}

int DemandCoveringProblem::period_to_execution_period(robust_model::period_id id) const {
    return int(id) / int(num_periods_per_planning_period());
}

void DemandCoveringProblem::for_each_planning_period_do(std::function<void(int)> f) const {
    for (size_t i = 0; i < num_planning_periods(); ++i) {
        f(i);
    }
}

void DemandCoveringProblem::for_each_execution_period_do(std::function<void(robust_model::period_id)> f) const {
    for_each_planning_period_do([&](size_t d) {
        for (robust_model::period_id i = 0; i < num_execution_periods_per_planning_period(); ++i) {
            f(robust_model::period_id(i + d * num_periods_per_planning_period() + 1));
        }
    });
}

void DemandCoveringProblem::for_each_location_do(std::function<void(Location::Index)> f) const {
    for (auto const location: location_ids()) {
        f(location);
    }
}

void DemandCoveringProblem::for_each_day_location_do(std::function<void(int, Location::Index)> f) const {
    for_each_planning_period_do([&](size_t d) {
        for_each_location_do([&](Location::Index l) {
            f(d, l);
        });
    });
}

void
DemandCoveringProblem::for_each_period_location_do(std::function<void(robust_model::period_id, Location::Index)> f) const {
    for_each_execution_period_do([&](robust_model::period_id i) {
        for_each_location_do([&](Location::Index l) {
            f(i, l);
        });
    });
}

void DemandCoveringProblem::for_each_period_location_pair_do(
        std::function<void(robust_model::period_id, Location::Index, Location::Index)> f) const {
    for_each_period_location_do([&](robust_model::period_id i, Location::Index l1) {
        for_each_location_do([&](Location::Index l2) {
            if (l1 == l2)
                return;
            f(i, l1, l2);
        });
    });
}

robust_model::period_id DemandCoveringProblem::previous_execution_period(robust_model::period_id period) const {
    return period - 1 - int((period - 1) % num_periods_per_planning_period() == 0);
}

robust_model::period_id DemandCoveringProblem::next_execution_period(robust_model::period_id period) const {
    return period + 1 + int((period + 1) % num_periods_per_planning_period() == 0);
}

double DemandCoveringProblem::lost_demand_cost() const {
    return _lost_demand_cost;
}

robust_model::ROModel& DemandCoveringROModelBuilder::ro_model() {
    return *_model_ptr;
}

DemandCoveringProblem const& DemandCoveringROModelBuilder::problem() {
    return _problem;
}

std::unique_ptr<robust_model::ROModel> DemandCoveringROModelBuilder::model_unique_ptr() {
    return std::move(_model_ptr);
}

DemandCoveringROModelBuilder::DemandCoveringROModelBuilder(DemandCoveringProblem const& problem)
        : _problem(problem),
          _model_ptr(std::make_unique<robust_model::ROModel>("Vaccination Problem")) {}

void DemandCoveringROModelBuilder::build() {
    add_uncertainty_variables();
    add_decision_variables();
    add_constraints();
    add_objective();
}

void DemandCoveringROModelBuilder::add_uncertainty_variables() {
    problem().for_each_period_location_do([&](robust_model::period_id period, Location::Index location) {
        _uncertainty_long.emplace(std::make_tuple(period, location),
                                  ro_model().add_uncertainty_variable("U_long_P" + std::to_string(period)
                                                                      + "_L" + std::to_string(location.raw_id()),
                                                                      problem().planning_period_period_id(
                                                                              problem().period_to_execution_period(
                                                                                      period)),
                                                                      0, 1));
        _uncertainty_short.emplace(std::make_tuple(period, location),
                                   ro_model().add_uncertainty_variable("U_short_P" + std::to_string(period)
                                                                       + "_L" + std::to_string(location.raw_id()),
                                                                       period,
                                                                       0, 1));
    });
}

void DemandCoveringROModelBuilder::add_decision_variables() {
    add_decision_variables_capacities();
    add_decision_variables_assignment();
    add_decision_variables_next_period_shift();
    add_decision_variables_other_location_shift();
}

void DemandCoveringROModelBuilder::add_decision_variables_capacities() {
    _capacity = ro_model().add_decision_variable("Capacity", -1, 0);
}

void DemandCoveringROModelBuilder::add_decision_variables_assignment() {
    problem().for_each_day_location_do([&](int day, Location::Index location) {
        _assignment.emplace(std::make_tuple(day, location),
                            ro_model().add_decision_variable("Assignment_D" + std::to_string(day)
                                                             + "_L" + std::to_string(location.raw_id()),
                                                             problem().planning_period_period_id(day),
                                                             0));
    });
}

void DemandCoveringROModelBuilder::add_decision_variables_next_period_shift() {
    problem().for_each_period_location_do([&](robust_model::period_id period, Location::Index location) {
        _next_period_shift.emplace(std::make_tuple(period, location),
                                   ro_model().add_decision_variable("Next_Period_Shift_P" + std::to_string(period)
                                                                    + "_L" + std::to_string(location.raw_id()),
                                                                    period,
                                                                    0));
    });
}

void DemandCoveringROModelBuilder::add_decision_variables_other_location_shift() {
    problem().for_each_period_location_pair_do(
            [&](robust_model::period_id period, Location::Index origin, Location::Index destination) {
                _other_location_shift.emplace(std::make_tuple(period, origin, destination),
                                              ro_model().add_decision_variable(
                                                      "OtherLocation_Shift_P" + std::to_string(period)
                                                      + "_O" + std::to_string(origin.raw_id())
                                                      + "_D" + std::to_string(destination.raw_id()),
                                                      period,
                                                      0));
            });
}

void DemandCoveringROModelBuilder::add_constraints() {
    add_constraints_capacity();
    add_constraints_demand_fulfillment();
}

void DemandCoveringROModelBuilder::add_constraints_capacity() {
    problem().for_each_planning_period_do([&](int day) {
        robust_model::RoAffineExpression sum;
        problem().for_each_location_do([&](Location::Index location) {
            sum += assignment(problem().planning_period_period_id(day), location);
        });
        ro_model().add_constraint(sum <= capacity(), "AssignmentCap_D" + std::to_string(day));
    });
}

void DemandCoveringROModelBuilder::add_constraints_demand_fulfillment() {
    problem().for_each_period_location_do([&](robust_model::period_id period, Location::Index location) {
        robust_model::RoAffineExpression lhs;
        lhs += assignment(period, location);
        lhs += next_period_shift(period, location);
        lhs -= (1 - problem().next_period_loss_factor(problem().previous_execution_period(period))) *
               next_period_shift(problem().previous_execution_period(period), location);
        problem().for_each_location_do([&](Location::Index other_location) {
            if (location == other_location)
                return;
            lhs += other_location_shift(period, location, other_location);
            lhs -= (1 - problem().distance_loss_factor(other_location->distance_to(location))) *
                   other_location_shift(period, other_location, location);
        });
        ro_model().add_constraint(lhs >= demand(period, location),
                                  "Demand_P" + std::to_string(period) + "_L" + std::to_string(location.raw_id()));
    });
}

void DemandCoveringROModelBuilder::add_objective() {
    robust_model::RoAffineExpression obj;
    obj += problem().capacity_cost() * capacity();
    problem().for_each_day_location_do([&](int day, Location::Index location) {
        obj += problem().capacity_usage_cost() * assignment(day, location);
    });
    problem().for_each_period_location_do([&](robust_model::period_id period, Location::Index location) {
        obj += problem().next_period_loss_factor(period) *
               next_period_shift(period, location) *
               problem().lost_demand_cost();
    });
    problem().for_each_period_location_pair_do(
            [&](robust_model::period_id period, Location::Index origin, Location::Index destination) {
                obj += problem().distance_loss_factor(origin->distance_to(destination)) *
                       problem().lost_demand_cost() *
                       other_location_shift(period, origin, destination);
            });
    ro_model().set_objective(obj, robust_model::ObjectiveSense::MIN);
}

robust_model::RoAffineExpression
DemandCoveringROModelBuilder::demand(robust_model::period_id period_id, Location::Index location) {
    return robust_model::RoAffineExpression(
            problem().base_demand(period_id, location) * (
                    1 +
                    _uncertainty_long.at(std::make_tuple(period_id, location)) +
                    .5 * _uncertainty_short.at(std::make_tuple(period_id, location)))
    );
}

robust_model::RoAffineExpression DemandCoveringROModelBuilder::capacity() {
    return robust_model::RoAffineExpression(_capacity.value());
}

robust_model::RoAffineExpression
DemandCoveringROModelBuilder::assignment(robust_model::period_id period, Location::Index location) {
    return robust_model::RoAffineExpression(
            _assignment.at(std::make_tuple(problem().period_to_execution_period(period), location)));
}

robust_model::RoAffineExpression
DemandCoveringROModelBuilder::next_period_shift(robust_model::period_id period, Location::Index location) {
    if (period == -1)
        return robust_model::RoAffineExpression(0);
    return robust_model::RoAffineExpression(_next_period_shift.at(std::make_tuple(period, location)));
}

robust_model::RoAffineExpression
DemandCoveringROModelBuilder::other_location_shift(robust_model::period_id period, Location::Index origin,
                                                   Location::Index destination) {
    return robust_model::RoAffineExpression(_other_location_shift.at(std::make_tuple(period, origin, destination)));
}

}
