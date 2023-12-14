//
// Created by simon on 25/10/2021.
//

#ifndef ROBUSTOPTIMIZATION_DEMANDCOVERINGINSTANCEGENERATOR_H
#define ROBUSTOPTIMIZATION_DEMANDCOVERINGINSTANCEGENERATOR_H

#include "../../test_helpers/InstanceGeneratorBase.h"

namespace testing {

class DemandCoveringInstanceGenerator : public InstanceGeneratorBase{
public:
    void add_num_locations(size_t m);

    void add_num_periods(size_t m);

    void add_lost_demand_cost(double c);

    std::string descriptions_test_specific() const override;

private:
    std::unique_ptr<robust_model::ROModel> generate_instance() override;

    std::string instance_description_test_specific() override;

    bool increment_test_specific() override;

private:
    std::vector<size_t> _numbers_of_locations;
    std::vector<size_t> _numbers_of_periods;
    std::vector<double> _lost_demand_costs;

    size_t _numbers_of_locations_id = 0,
            _numbers_of_periods_id = 0,
            _lost_demand_costs_id = 0;

};

}

#endif //ROBUSTOPTIMIZATION_DEMANDCOVERINGINSTANCEGENERATOR_H
