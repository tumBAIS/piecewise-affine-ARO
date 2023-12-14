//
// Created by simon on 12.05.21.
//

#include <random>
#include <utility>
#include <set>
#include "UncertaintySet.h"
#include "ROModel.h"
#include "SOCModel.h"
#include "../helpers/helpers.h"

namespace robust_model {

UncertaintyRealization::UncertaintyRealization(std::vector<double> values) : _values(std::move(values)) {}

std::vector<double> const& UncertaintyRealization::values() const {
    return _values;
}

double UncertaintyRealization::value(UncertaintyVariable::Index const& id) const {
    return _values.at(id.raw_id());
}

std::string UncertaintySet::to_string(UncertaintySet::SpecialSetType set_type) {
    switch (set_type) {
        case SpecialSetType::BOX:
            return "BOX";
        case SpecialSetType::BALL:
            return "BALL";
        case SpecialSetType::BUDGET:
            return "BUDGET";
        case SpecialSetType::OTHER:
            return "OTHER";
    }
    helpers::exception_check(false, "Illegal set special_type!");
    return "";
}

UncertaintySetConstraintsSet::UncertaintySetConstraintsSet(
        UncertaintySetConstraintsSet::Index id, UncertaintySet const& uncertainty_set) :
        IndexedObject(id), _uncertainty_set(uncertainty_set) {}

void
UncertaintySetConstraintsSet::add_uncertainty_constraint(UncertaintySetConstraintsSet::Constraint const& constraint) {
    _uncertainty_constraints.emplace_back(constraint);
}

std::vector<UncertaintySetConstraintsSet::Constraint> const& UncertaintySetConstraintsSet::constraints() const {
    return _uncertainty_constraints;
}

bool UncertaintySetConstraintsSet::is_box() const {
    return std::all_of(_uncertainty_constraints.begin(), _uncertainty_constraints.end(),
                       [](auto const& constr) {
                           return constr.is_variable_bound();
                       });
}

void UncertaintySetConstraintsSet::clear() {
    _uncertainty_constraints.clear();
}

bool UncertaintySetConstraintsSet::is_empty() const {
    return constraints().empty();
}

bool UncertaintySetConstraintsSet::is_norm_ball() const {
    if (constraints().size() != 1)
        return false;
    auto const& constr = constraints().front();
    if (constr.sense() != ConstraintSense::LEQ)
        return false;
    if (not constr.expression().affine().linear().empty())
        return false;
    if (constr.expression().affine().constant() > 0)
        return false;
    if (not constr.expression().normed_vector().is_variable_vector_norm(_uncertainty_set.indices()))
        return false;
    return true;
}

double UncertaintySetConstraintsSet::budget() const {
    helpers::exception_check(is_norm_ball(), "Can only find budget for norm balls!");
    return - constraints().front().expression().affine().constant();
}


UncertaintySet::UncertaintySet(ROModel const& model) : _model(model) {
    helpers::IndexedObjectOwner<UncertaintySetConstraintsSet>::base_add_object(*this);
}

UncertaintySet::UncertaintyReference
UncertaintySet::add_variable(std::string const& name, std::optional<period_id> p, double lb, double ub) {
    return UncertaintyReference(helpers::IndexedObjectOwner<UncertaintyVariable>::base_add_object(name, p, lb, ub));
}

std::vector<UncertaintyVariable> const&
UncertaintySet::variables() const {
    return helpers::IndexedObjectOwner<UncertaintyVariable>::objects();
}

std::vector<double> UncertaintySet::lower_bounds() const {
    std::vector<double> lbs(num_variables());
    for (auto const& var: variables()) {
        lbs[var.id().raw_id()] = var.lb();
    }
    return lbs;
}

std::vector<double> UncertaintySet::upper_bounds() const {
    std::vector<double> ubs(num_variables());
    for (auto const& var: variables()) {
        ubs[var.id().raw_id()] = var.ub();
    }
    return ubs;
}

size_t UncertaintySet::num_variables() const {
    return helpers::IndexedObjectOwner<UncertaintyVariable>::num_objects();
}

ROModel const&
robust_model::UncertaintySet::model() const {
    return _model;
}

UncertaintySet::SpecialSetType
robust_model::UncertaintySet::special_type() const {
    if (constraint_sets().size() > 1)
        return SpecialSetType::OTHER;
    auto const& constr_set = constraint_sets().front();
    if (constr_set->is_empty() and all_bounded_variables())
        return SpecialSetType::BOX;
    if (constr_set->is_norm_ball()) {
        if (constr_set->constraints().front().expression().normed_vector().norm_type() == VectorNormType::Two)
            return SpecialSetType::BALL;
        if (constr_set->constraints().front().expression().normed_vector().norm_type() == VectorNormType::One)
            return SpecialSetType::BUDGET;
    }
    return UncertaintySet::SpecialSetType::OTHER;
}

void
UncertaintySet::add_special_type_constraint(UncertaintySet::SpecialSetType type, double budget) {
    switch (type) {
        case SpecialSetType::BALL:
            add_uncertainty_constraint({SOCExpression<UncertaintyVariable>::norm(
                    variables(), VectorNormType::Two) <= budget,
                                        "BallConstraint"});
            return;
        case SpecialSetType::BUDGET:
            add_uncertainty_constraint({SOCExpression<UncertaintyVariable>::norm(
                    variables(), VectorNormType::One) <= budget,
                                        "BudgetConstraint"});
            return;
        default:
            helpers::exception_throw("No special Constraint or not implemented!");
    }
}

double
UncertaintySet::budget() const {
    helpers::exception_check(constraint_sets().size() == 1,
                             "Can only find budget for single norm balls!"
                             );
    return constraint_sets().front()->budget();
}

void UncertaintySet::add_uncertainty_constraint(UncertaintySet::Constraint const& constraint) {
    add_uncertainty_constraint(constraint, constraint_sets().back());
}

void UncertaintySet::add_uncertainty_constraint(UncertaintySet::Constraint const& constraint,
                                                UncertaintySetConstraintsSet::Index constraint_set) {
    helpers::IndexedObjectOwner<UncertaintySetConstraintsSet>::object(constraint_set).add_uncertainty_constraint(
            constraint);
}

std::vector<UncertaintySet::Constraint> const&
UncertaintySet::uncertainty_constraints() const {
    helpers::exception_check(constraint_sets().size() == 1,
                             "Only use the constraints funtion, when there is just one set of constraints");
    return uncertainty_constraints(helpers::IndexedObjectOwner<UncertaintySetConstraintsSet>::ids().front());
}

std::vector<UncertaintySet::Constraint> const&
UncertaintySet::uncertainty_constraints(UncertaintySetConstraintsSet::Index constraint_set) const {
    return constraint_set->constraints();
}

UncertaintySetConstraintsSet::Index UncertaintySet::add_constraint_set() {
    return helpers::IndexedObjectOwner<UncertaintySetConstraintsSet>::base_add_object(*this);
}

std::vector<UncertaintySetConstraintsSet::Index> const& UncertaintySet::constraint_sets() const {
    return helpers::IndexedObjectOwner<UncertaintySetConstraintsSet>::ids();
}

void UncertaintySet::clear_constraints() {
    helpers::IndexedObjectOwner<robust_model::UncertaintySetConstraintsSet>::clear();
    add_constraint_set();
}

std::string
UncertaintySet::full_string() const {
    std::string s = "Uncertainty Set of " + model().name() + "\n";
    s += to_string(special_type());
    if (special_type() == SpecialSetType::OTHER) {
        for (auto const& constr: uncertainty_constraints()) {
            s += "\n" + constr.to_string();
        }
    }
    for (auto const& var: helpers::IndexedObjectOwner<UncertaintyVariable>::objects()) {
        s += "\n" + std::to_string(var.lb()) + " <= " + var.name() + " p(" + std::to_string(var.period()) + ") <= " +
             std::to_string(var.ub());
    }
    return s;
}

std::vector<std::vector<double>> UncertaintySet::generate_uncertainty(size_t const num_realizations) const {
    switch (special_type()) {
        case SpecialSetType::BALL:
            return generate_uncertainty_ball(num_realizations);
        case SpecialSetType::BUDGET:
            return generate_uncertainty_budgeted(num_realizations);
        default:
            helpers::exception_check(false, "Only BALL and BUDGET implemented yet!");
    }
    return {};
}

std::vector<std::vector<double>> UncertaintySet::generate_uncertainty_ball(size_t num_realizations) const {
    helpers::exception_check(all_0_1_variables(),
                             "Generation of uncertain samples only implemented for [0,1] uncertainty variables yet!");
    helpers::exception_check(budget() == 1, "BALL uncertainty generation only implemented for radius 1!");
    std::ranlux48 generator(std::chrono::system_clock::now().time_since_epoch().count());
    std::normal_distribution<double> normal_distribution(0., 1.);
    std::uniform_real_distribution<double> uniform_distribution(0., 1.);

    std::vector<std::vector<double>> realizations(num_realizations, std::vector<double>(num_variables()));
    for (size_t i = 0; i < num_realizations; ++i) {
        auto& realization = realizations.at(i);
        double norm = 0;
        for (size_t j = 0; j < num_variables(); ++j) {
            double rn = std::abs(normal_distribution(generator));
            realization[j] = rn;
            norm += rn * rn;
        }
        norm = std::sqrt(norm);
        double const radius = std::pow(uniform_distribution(generator), 1. / double(num_variables()));
        double const scale = radius / norm;
        for (size_t j = 0; j < num_variables(); ++j) {
            realization[j] *= scale;
        }
    }
    return realizations;
}

std::vector<std::vector<double>> UncertaintySet::generate_uncertainty_budgeted(size_t num_realizations) const {
    helpers::exception_check(all_0_1_variables(),
                             "Generation of uncertain samples only implemented for [0,1] uncertainty variables yet!");
    std::ranlux48 generator(std::chrono::system_clock::now().time_since_epoch().count());
    std::exponential_distribution<double> exponential_distribution(1.);

    std::vector<std::vector<double>> realizations(num_realizations, std::vector<double>(num_variables()));
    bool previous_in_set = true;
    for (size_t i = 0; i < num_realizations; i += previous_in_set) {
        previous_in_set = true;
        auto& realization = realizations.at(i);
        double one_norm = exponential_distribution(generator);
        for (size_t j = 0; j < num_variables(); ++j) {
            double rn = exponential_distribution(generator);
            realization[j] = rn;
            one_norm += rn;
        }
        double const scale = budget() / one_norm;
        for (size_t j = 0; j < num_variables(); ++j) {
            realization[j] *= scale;
            if (realization[j] > 1) {
                previous_in_set = false;
                break;
            }
        }
    }
    helpers::exception_check(previous_in_set, "When done, the last realization should be valid!");
    return realizations;
}

std::vector<UncertaintyVariable::Index> const& UncertaintySet::indices() const {
    return helpers::IndexedObjectOwner<UncertaintyVariable>::ids();
}

std::unique_ptr<SOCModel> UncertaintySet::to_soc_model() const {
    auto soc_model_ptr = std::make_unique<SOCModel>("Uncertainty" + model().name());
    auto& soc_model = *soc_model_ptr;
    std::vector<SOCVariable::Reference> vars;
    for (auto const& var: variables()) {
        vars.emplace_back(soc_model.add_variable(var.name(), var.lb(), var.ub()));
    }
    for (auto const& constr: uncertainty_constraints()) {
        soc_model.add_constraint(robust_model::SOCModel::Constraint(
                constr.sense(),
                constr.expression().translate_to_other(vars),
                constr.name()));
    }
    return soc_model_ptr;
}

bool UncertaintySet::all_0_1_variables() const {
    return std::all_of(variables().begin(), variables().end(),
                       [](UncertaintyVariable const& v) { return v.is_0_1_var(); });
}

bool UncertaintySet::all_bounded_variables() const {
    return std::all_of(variables().begin(), variables().end(),
                       [](UncertaintyVariable const& v) { return v.bounded(); });
}

bool UncertaintySet::rotational_invariant() const {
    if (not std::set{SpecialSetType::BOX,
                     SpecialSetType::BUDGET,
                     SpecialSetType::BALL}.contains(special_type())) {
        helpers::warning_throw("Rotational Invariance is only checked for special sets!");
        return false;
    }
    if (variables().empty()) {
        return true;
    }
    double const lb = variables().front().lb();
    double const ub = variables().front().ub();
    return std::all_of(variables().begin(), variables().end(),
                       [lb, ub](auto const& var) { return var.lb() == lb and var.ub() == ub; });
}

bool UncertaintySet::symmetric() const {
    if (not std::set{SpecialSetType::BOX,
                     SpecialSetType::BUDGET,
                     SpecialSetType::BALL}.contains(special_type())) {
        helpers::warning_throw("Symmetry is only checked for special sets!");
        return false;
    }
    if (variables().empty()) {
        return true;
    }
    return std::all_of(variables().begin(), variables().end(),
                       [](auto const& var) { return -var.lb() == var.ub(); });
}

bool UncertaintySet::non_negative() const {
    return std::all_of(variables().begin(), variables().end(),
                       [](auto const& var) { return var.lb() >= 0; });
}

double UncertaintySet::max_one_norm_k_active(size_t k) const {
    if (k == 0) {
        return 0;
    }
    helpers::exception_check(k <= num_variables(), "Can't have more then num_var active variables!");
    if (rotational_invariant()) {
        double max_bound = std::max(-variables().front().lb(), variables().front().ub());
        double const box_bound = double(k) * max_bound;

        switch (special_type()) {
            case SpecialSetType::BOX:
                return box_bound;
            case SpecialSetType::BUDGET:
                return std::min(box_bound, budget());
            case SpecialSetType::BALL:
                return std::min(box_bound, std::sqrt(double(k)));
            default:
                helpers::exception_throw("Not implemented!");
        }
    }
    helpers::exception_throw("Not implemented!");
}

}
