//
// Created by simon on 20/10/2021.
//

#ifndef ROBUSTOPTIMIZATION_DEMANDCOVERINGPROBLEM_H
#define ROBUSTOPTIMIZATION_DEMANDCOVERINGPROBLEM_H

#include "../../../helpers/helpers.h"
#include "../../../solvers/aro_policy_solvers/AffineAdjustablePolicySolver.h"
#include <functional>
#include <memory>

namespace demand_covering_problem {

class Location: public helpers::IndexedObject<Location>{
public:
    Location(helpers::SmartIndex<Location> const& id, double x_coord, double y_coord);

    double distance_to(Location::Index other) const;

    double x_coord() const;

    double y_coord() const;

private:
    double _x_coord;
    double _y_coord;
};


class DemandCoveringProblem: helpers::IndexedObjectOwner<Location> {
public:
    DemandCoveringProblem(size_t num_planning_periods, size_t num_execution_per_planning_periods,
                          double capacity_cost, double capacity_usage_cost, double distance_loss_factor,
                          double next_period_loss, double next_day_loss, double lost_demand_cost);

    helpers::SmartIndex<Location> add_location(double x, double y);

    void set_demand_pattern(std::function<double(robust_model::period_id, Location::Index)> demand_pattern);

    double base_demand(robust_model::period_id period, Location::Index location) const;

    size_t num_locations() const;

    std::vector<Location> const& locations() const;

    std::vector<Location::Index> const& location_ids() const;

    size_t num_planning_periods() const;

    size_t num_execution_periods_per_planning_period() const;

    size_t num_periods_per_planning_period() const;

    double capacity_cost() const;

    double capacity_usage_cost() const;

    double lost_demand_cost() const;

    double distance_loss_factor(double dist) const;

    double next_period_loss_factor(robust_model::period_id period) const;

    robust_model::period_id planning_period_period_id(size_t day) const;

    int period_to_execution_period(robust_model::period_id id) const;

    robust_model::period_id previous_execution_period(robust_model::period_id period) const;
    robust_model::period_id next_execution_period(robust_model::period_id period) const;

    void for_each_planning_period_do(std::function<void(int)> f) const;
    void for_each_execution_period_do(std::function<void(robust_model::period_id)> f) const;
    void for_each_location_do(std::function<void(Location::Index)> f) const;
    void for_each_day_location_do(std::function<void(int, Location::Index)> f) const;
    void for_each_period_location_do(std::function<void(robust_model::period_id , Location::Index)> f) const;
    void for_each_period_location_pair_do(std::function<void(robust_model::period_id, Location::Index , Location::Index)> f) const;

private:
    size_t const _num_planning_periods;
    size_t const _num_execution_per_planning_periods;

    double const _capacity_cost;
    double const _capacity_usage_cost;

    double const _distance_loss_factor;
    double const _next_period_loss;
    double const _next_day_loss;

    double const _lost_demand_cost;

    std::function<double(robust_model::period_id, Location::Index)> _base_demand_pattern;
};

class DemandCoveringROModelBuilder{
public:
    explicit DemandCoveringROModelBuilder(DemandCoveringProblem const& problem);

    void build();

    std::unique_ptr<robust_model::ROModel> model_unique_ptr();
private:
    robust_model::ROModel& ro_model();
    DemandCoveringProblem const& problem();

    void add_uncertainty_variables();

    void add_decision_variables();
    void add_decision_variables_capacities();
    void add_decision_variables_assignment();
    void add_decision_variables_next_period_shift();
    void add_decision_variables_other_location_shift();

    void add_constraints();
    void add_constraints_capacity();
    void add_constraints_demand_fulfillment();

    void add_objective();

    robust_model::RoAffineExpression demand(robust_model::period_id period_id, Location::Index location);

    robust_model::RoAffineExpression capacity();
    robust_model::RoAffineExpression assignment(robust_model::period_id period, Location::Index location);
    robust_model::RoAffineExpression next_period_shift(robust_model::period_id period, Location::Index location);
    robust_model::RoAffineExpression other_location_shift(robust_model::period_id period, Location::Index origin, Location::Index destination);
private:
    std::map<std::tuple<robust_model::period_id , Location::Index>, robust_model::UncertaintyVariable::Reference> _uncertainty_long;
    std::map<std::tuple<robust_model::period_id , Location::Index>, robust_model::UncertaintyVariable::Reference> _uncertainty_short;

    std::optional<robust_model::DecisionVariable::Reference> _capacity;
    std::map<std::tuple<int, Location::Index>, robust_model::DecisionVariable::Reference> _assignment;
    std::map<std::tuple<robust_model::period_id , Location::Index>, robust_model::DecisionVariable::Reference> _next_period_shift;
    std::map<std::tuple<robust_model::period_id, Location::Index, Location::Index>, robust_model::DecisionVariable::Reference> _other_location_shift;


private:
    std::unique_ptr<robust_model::ROModel> _model_ptr;
    DemandCoveringProblem const& _problem;
};

}

#endif //ROBUSTOPTIMIZATION_DEMANDCOVERINGPROBLEM_H
