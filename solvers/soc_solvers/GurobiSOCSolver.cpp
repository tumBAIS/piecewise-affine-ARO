//
// Created by thomae on 12/8/2023.
//

#include "GurobiSOCSolver.h"
#include "../../models/SOCModel.h"

namespace solvers {

GurobiSOCSolver::GurobiSOCSolver(robust_model::SOCModel& soc_model) :
        SOCSolverBase(soc_model) {}


void GurobiSOCSolver::solve_implementation() {
    try {
        set_status(SolverBase::Status::UNSOLVED);
        _grb_model->optimize();
        if (gurobi_model().get(GRB_IntAttr_Status) == GRB_OPTIMAL)
            set_status(SolverBase::Status::OPTIMAL);
        if (gurobi_model().get(GRB_IntAttr_Status) == GRB_TIME_LIMIT)
            set_status(SolverBase::Status::TIME_LIMIT);
    } catch (GRBException e) {
        if (e.getErrorCode() == GRB_ERROR_OUT_OF_MEMORY) {
            set_status(SolverBase::Status::MEMORY_LIMIT);
        } else {
            helpers::exception_check(false, "Gurobi error "
                                            + std::to_string(e.getErrorCode()) + ": " + e.getMessage());
        }
    }
    if (has_solution()) {
        set_runtime(gurobi_model().get(GRB_DoubleAttr_Runtime));
        set_objective_value(gurobi_model().get(GRB_DoubleAttr_ObjVal));
    }
    helpers::warning_check(
            status() != SolverBase::Status::UNSOLVED,
            "No gurobi solution found! Unexpected Return Code " +
            std::to_string(gurobi_model().get(GRB_IntAttr_Status)));

    transfer_solution_to_soc_model();
}

void GurobiSOCSolver::build_implementation() {
    _grb_env = std::make_unique<GRBEnv>(true);
    _grb_env->set(GRB_IntParam_LogToConsole, 0);
    _grb_env->set(GRB_IntParam_Threads, 1);
    _grb_env->set(GRB_DoubleParam_MIPGap, 0.01);
    _grb_env->start();
    _grb_model = std::make_unique<GRBModel>(*_grb_env);
    _grb_vars = std::make_unique<std::vector<GRBVar>>();
}

void GurobiSOCSolver::update_implementation() {
    if (has_runtime_limit())
        gurobi_model().set(GRB_DoubleParam_TimeLimit, runtime_limit());
    if (has_memory_limit())
        gurobi_model().set(GRB_DoubleParam_MemLimit, memory_limit());
    update_variables();
    update_soc_constraints();
    update_sos_constraints();
    update_objectives();

    gurobi_model().update();
}

void GurobiSOCSolver::update_variables() {
    for (auto const& var: variables_to_add_grb()) {
        _grb_vars->emplace_back(gurobi_model().addVar(var.lb(), var.ub(), 0, to_grb_type(var.type())));
    }
    _grb_next_var_to_add = soc_model().variables().size();
    gurobi_model().update();
}

void GurobiSOCSolver::update_soc_constraints() {
    for (auto const& constr: soc_constraints_to_add_grb()) {
        if (constr.soc_expression().is_affine()) {
            gurobi_model().addConstr(to_gurobi_linear(constr.soc_expression().affine()),
                                     to_grb_sense(constr.sense()), 0,
                                     constr.name()
            );
            continue;
        }
        switch (constr.soc_expression().normed_vector().norm_type()) {
            case robust_model::VectorNormType::Two: {
                auto expr = GRBQuadExpr();
                for (auto const& affine: constr.soc_expression().normed_vector().normed_vector()) {
                    auto const lin_expr = to_gurobi_linear(affine);
                    expr += lin_expr * lin_expr;
                }
                auto const lin_expr = to_gurobi_linear(constr.soc_expression().affine());
                gurobi_model().addQConstr(expr, to_grb_sense(constr.sense()), lin_expr * lin_expr, constr.name());
                gurobi_model().addConstr(lin_expr <= 0, constr.name() + "_pos");
                continue;
            }
            case robust_model::VectorNormType::One: {
                auto expr = GRBLinExpr();
                for (auto const& affine: constr.soc_expression().normed_vector().normed_vector()) {
                    auto const lin_expr = to_gurobi_linear(affine);
                    auto const abs = gurobi_model().addVar(0, robust_model::NO_VARIABLE_UB, 0, GRB_CONTINUOUS);
                    gurobi_model().addConstr(abs >= lin_expr);
                    gurobi_model().addConstr(abs >= -lin_expr);
                    expr += abs;
                }
                auto const lin_expr = to_gurobi_linear(constr.soc_expression().affine());
                gurobi_model().addConstr(expr, to_grb_sense(constr.sense()), -lin_expr, constr.name());
                continue;
            }
            case robust_model::VectorNormType::Max: {
                auto const lin_expr = to_gurobi_linear(constr.soc_expression().affine());
                for (auto const& affine: constr.soc_expression().normed_vector().normed_vector()) {
                    gurobi_model().addConstr(to_gurobi_linear(affine), to_grb_sense(constr.sense()), -lin_expr);
                    gurobi_model().addConstr(-to_gurobi_linear(affine), to_grb_sense(constr.sense()), -lin_expr);
                }
                continue;
            }
        }
    }
    _grb_next_constr_to_add = soc_model().soc_constraints().size();
}

void GurobiSOCSolver::update_sos_constraints() {
    for (auto const& constr: sos_constraints_to_add_grb()) {
        size_t const len = constr.exclusive_variables().size();
        std::vector<GRBVar> exclusive_vars;
        std::vector<double> weights;
        for (auto const& var: constr.exclusive_variables()) {
            exclusive_vars.emplace_back(_grb_vars->at(var.raw_id()));
            weights.emplace_back(double(weights.size()));
        }
        gurobi_model().addSOS(exclusive_vars.data(), weights.data(), int(len), GRB_SOS_TYPE1);
    }
    _grb_next_sos_constr_to_add = soc_model().sos_constraints().size();
}

void GurobiSOCSolver::update_objectives() {
    static const int max_objectives = 100;
    if (soc_model().is_multi_objective()) {
        gurobi_model().set(GRB_IntAttr_ModelSense, to_grb_sense(soc_model().objective().sense()));
        gurobi_model().set(GRB_IntAttr_NumObj, int(soc_model().objectives().size()));
        for (size_t i = _grb_next_obj_to_add; i < soc_model().objectives().size(); ++i) {
            helpers::exception_check(
                    i < max_objectives, "Did not expect to get more than 100 objectives! "
                                        "Adjust max_objectives when more are needed!");
            gurobi_model().setObjectiveN(to_gurobi_linear(soc_model().objective(i).expression()), int(i),
                                         int(max_objectives - i));
        }
    } else if (_grb_next_obj_to_add == 0) {
        gurobi_model().setObjective(to_gurobi_linear(soc_model().objective().expression()),
                                    to_grb_sense(soc_model().objective().sense()));
    }
    _grb_next_obj_to_add = soc_model().objectives().size();
}

GRBLinExpr GurobiSOCSolver::to_gurobi_linear(
        robust_model::AffineExpression<robust_model::SOCVariable::Reference> const& affine) const {
    GRBLinExpr expr(affine.constant());
    for (auto const& svar: affine.linear().scaled_variables()) {
        expr += svar.scale() * _grb_vars->at(svar.variable().raw_id());
    }
    return expr;;
}

helpers::VectorSlice<robust_model::SOCVariable> GurobiSOCSolver::variables_to_add_grb() const {
    return {soc_model().variables(), _grb_next_var_to_add, soc_model().variables().size()};
}

helpers::VectorSlice<robust_model::SOCConstraint<robust_model::SOCVariable>>
GurobiSOCSolver::soc_constraints_to_add_grb() const {
    return {soc_model().soc_constraints(), _grb_next_constr_to_add, soc_model().soc_constraints().size()};
}

helpers::VectorSlice<robust_model::SOSConstraint> GurobiSOCSolver::sos_constraints_to_add_grb() const {
    return {soc_model().sos_constraints(), _grb_next_sos_constr_to_add, soc_model().sos_constraints().size()};
}

helpers::VectorSlice<robust_model::SOCModel::Objective> GurobiSOCSolver::objectives_to_add_grb() const {
    return {soc_model().objectives(), _grb_next_obj_to_add, soc_model().objectives().size()};
}

GRBModel& GurobiSOCSolver::gurobi_model() {
    return *_grb_model;
}

void GurobiSOCSolver::transfer_solution_to_soc_model() {
    if (has_solution()) {
        non_const_soc_model().set_solution(*this);
        if (soc_model().all_affine() and (not soc_model().is_multi_objective()) and soc_model().is_continuous()) {
            std::vector<double> dual_values(soc_model().soc_constraints().size());
            for (int i = 0; i < dual_values.size(); ++i) {
                dual_values.at(i) = gurobi_model().getConstr(i).get(GRB_DoubleAttr_Pi);
            }
            non_const_soc_model().set_dual_soc_constraint_values(dual_values);
        }
    } else {
        non_const_soc_model().invalidate_solution();
    }
}

double GurobiSOCSolver::value(robust_model::SOCVariable::Index const& id) const {
    helpers::exception_check(has_solution(), "Only extract solution, when it exists!");
    return _grb_vars->at(id.raw_id()).get(GRB_DoubleAttr_X);
}

void GurobiSOCSolver::objectives_reset() {
    _grb_next_obj_to_add = 0;
}

}