//
// Created by simon on 20/05/2021.
//

#ifndef ROBUSTOPTIMIZATION_SOCMODEL_H
#define ROBUSTOPTIMIZATION_SOCMODEL_H

#include "basic_model_objects/SOCConstraint.h"
#include "basic_model_objects/SOCVariable.h"
#include "basic_model_objects/ObjectiveBase.h"
#include "../helpers/helpers.h"
#include "basic_model_objects/SOSConstraint.h"
#include <memory>

class GRBModel;
class GRBEnv;

namespace robust_model {

class SOCModel : public helpers::IndexedObjectOwner<SOCVariable> {
public:
    using SOCVariableReference = VariableReference<SOCVariable>;
    using Constraint = SOCConstraint<SOCVariable>;
    using Objective = ObjectiveBase<AffineExpression<SOCVariable::Reference>>;
public:
    explicit SOCModel(std::string  name);

    SOCVariableReference add_variable(std::string const& name, double lb = NO_VARIABLE_LB, double ub=NO_VARIABLE_UB, VariableType type = VariableType::Continuous);
    SOCVariableReference add_variable(std::string const& name, double lb, double ub, DualInformation const& info);

    std::vector<SOCVariable::Reference> add_variables(std::size_t n, std::string const& name, double lb = NO_VARIABLE_LB, double ub=NO_VARIABLE_UB, VariableType type = VariableType::Continuous);

    void add_constraint(SOCConstraint<SOCVariable> const& constraint);

    template<ConvertibleTo<SOCConstraint<SOCVariable>::RawConstraint> RC>
    void add_constraint(RC const& raw_constraint, std::string const& name){
        add_constraint({raw_constraint, name});
    }

    void add_sos_constraint(std::vector<SOCVariable::Reference> const& exclusive_variables);

    std::vector<SOCVariable> const& variables() const;

    std::vector<SOCConstraint<SOCVariable>> const& soc_constraints() const;

    std::vector<SOSConstraint> const& sos_constraints() const;

    Objective const& objective(size_t id = 0) const;

    std::vector<Objective> const& objectives() const;

    void clear_and_set_objective(Objective const& objective);

    void add_objective(Objective const& objective);

    bool is_multi_objective() const;

    bool is_continuous() const;

    bool all_affine() const;

    std::string const& name() const;

    std::string full_string() const;

    template<class S>
    void set_solution(S const& solution){
        for(auto & var : non_const_objects()){
            var.set_solution(solution.value(var.id()));
        }
        for (auto& objective: _objectives) {
            objective.set_solution(objective.expression().value(solution));
        }
    }

    void set_dual_soc_constraint_values(std::vector<double> const& values);

    void invalidate_solution();
    
    void compute_dual();
    SOCModel const& dual() const;
    SOCModel & dual();

private:
    void invalidate_dual();

    void
    add_dual_of_var_bound(AffineExpression<SOCVariable::Reference> & dual_objective,
                          AffineExpression<SOCVariable::Reference> & dual_constraint_expression,
                          SOCVariable::Index const& variable
                          );

    void
    add_dual_of_expression(AffineExpression<SOCVariable::Reference> & dual_objective,
                           std::vector<AffineExpression < SOCVariable::Reference>> & dual_constraint_expressions,
                           size_t constraint_number);

    void
    add_dual_of_normed_vector(AffineExpression<SOCVariable::Reference> & dual_objective,
                              std::vector<AffineExpression < SOCVariable::Reference>> & dual_constraint_expressions,
                              SOCVariable::Reference dv,
                              size_t constraint_number);

    void
    add_dual_of_objective(AffineExpression<SOCVariable::Reference> & dual_objective,
                           std::vector<AffineExpression <SOCVariable::Reference>> & dual_constraint_expressions);

private:
    std::string const _name;
    std::vector<Objective> _objectives;
    std::vector<SOCConstraint<SOCVariable>> _soc_constraints;
    std::vector<SOSConstraint> _sos_constraints;
    std::unique_ptr<SOCModel> _dual;
};

}
#endif //ROBUSTOPTIMIZATION_SOCMODEL_H
