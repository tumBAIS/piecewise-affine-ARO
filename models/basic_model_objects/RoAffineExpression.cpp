//
// Created by simon on 08.05.21.
//

#include "RoAffineExpression.h"
#include "SolutionRealization.h"

#include <utility>

namespace robust_model {

RoAffineExpression::RoAffineExpression(double constant, LinearExpression<DecisionVariable::Reference> const& decisions,
                                       LinearExpression<UncertaintyVariable::Reference> const& uncertainties,
                                       LinearExpression<UncertaintyScaledDecision> const& uncertainty_decisions,
                                       MultiUncertaintyBehaviour multi_uncertainty_behaviour)
        : _constant(constant),
          _decisions(decisions),
          _uncertainties(uncertainties),
          _uncertainty_decisions(uncertainty_decisions),
          _multi_uncertainty_behaviour(multi_uncertainty_behaviour){}

RoAffineExpression::RoAffineExpression(double constant, LinearExpression<DecisionVariable::Reference> const& decisions,
                                       LinearExpression<UncertaintyVariable::Reference> const& uncertainties,
                                       LinearExpression<UncertaintyScaledDecision> const& uncertainty_decisions)
        : _constant(constant),
          _decisions(decisions),
          _uncertainties(uncertainties),
          _uncertainty_decisions(uncertainty_decisions) {}

RoAffineExpression::RoAffineExpression(double constant,
                                       LinearExpression<DecisionVariable::Reference> const& decisions,
                                       LinearExpression<UncertaintyVariable::Reference> const& uncertainties)
        : _constant(
        constant), _decisions(decisions), _uncertainties(uncertainties) {}

RoAffineExpression::RoAffineExpression(double constant,
                                       std::vector<ScaledVariable<DecisionVariable::Reference>> const& decisions,
                                       std::vector<ScaledVariable<UncertaintyVariable::Reference>> const& uncertainties)
        : _constant(constant), _decisions{decisions}, _uncertainties{uncertainties} {}

RoAffineExpression::RoAffineExpression(AffineExpression<DecisionVariable::Reference> const& t) :
        _constant(t.constant()), _decisions(t.linear()) {}

RoAffineExpression::RoAffineExpression(AffineExpression<UncertaintyVariable::Reference> const& t) :
        _constant(t.constant()), _uncertainties(t.linear()) {}

RoAffineExpression::RoAffineExpression(AffineExpression<UncertaintyScaledDecision> const& t) :
        _constant(t.constant()), _uncertainty_decisions(t.linear()) {}

RoAffineExpression::RoAffineExpression(double c) : _constant(c) {}

RoAffineExpression& RoAffineExpression::operator*=(double scale) {
    _constant *= scale;
    _uncertainties *= scale;
    _decisions *= scale;
    _uncertainty_decisions *= scale;
    return *this;
}

RoAffineExpression& RoAffineExpression::operator/=(double scale) {
    return *this *= (1. / scale);
}

RoAffineExpression RoAffineExpression::operator*(double scale) const {
    return RoAffineExpression(*this) *= scale;
}

RoAffineExpression RoAffineExpression::operator/(double scale) const {
    return RoAffineExpression(*this) /= scale;
}

RoAffineExpression& RoAffineExpression::operator+=(RoAffineExpression const& other) {
    _constant += other.constant();
    _decisions += other.decisions();
    _uncertainties += other.uncertainties();
    _uncertainty_decisions += other.uncertainty_decisions();
    return *this;
}

RoAffineExpression& RoAffineExpression::operator+=(double c) {
    _constant += c;
    return *this;
}

RoAffineExpression& RoAffineExpression::operator+=(AffineExpression<UncertaintyScaledDecision> const& expr) {
    _constant += expr.constant();
    _uncertainty_decisions += expr.linear();
    return *this;
}


RoAffineExpression RoAffineExpression::operator-() const {
    return {-constant(), -decisions(),
            -uncertainties(), -uncertainty_decisions(),
            multi_uncertainty_behaviour()};
}

double RoAffineExpression::constant() const {
    return _constant;
}

LinearExpression<DecisionVariable::Reference> const& RoAffineExpression::decisions() const {
    return _decisions;
}

LinearExpression<UncertaintyVariable::Reference> const& RoAffineExpression::uncertainties() const {
    return _uncertainties;
}

LinearExpression<UncertaintyScaledDecision> const& RoAffineExpression::uncertainty_decisions() const {
    return _uncertainty_decisions;
}

std::string RoAffineExpression::to_string() const {
    return std::to_string(constant())
           + (decisions().empty() ? "" : " + ")
           + decisions().to_string()
           + (uncertainty_decisions().empty() ? "" : " + ")
           + uncertainty_decisions().to_string()
           + (uncertainties().empty() ? "" : " + ")
           + uncertainties().to_string();
}

double RoAffineExpression::value(SolutionRealization const& realization) const {
    return decisions().value(realization) + uncertainties().value(realization) + constant();
}

double RoAffineExpression::lb() const {
    return uncertainties().lb() + decisions().lb() + uncertainty_decisions().lb() + constant();
}

double RoAffineExpression::ub() const {
    return uncertainties().ub() + decisions().ub() + uncertainty_decisions().ub() + constant();
}

bool RoAffineExpression::is_zero() const {
    return uncertainties().empty() and
           decisions().empty() and
           uncertainty_decisions().empty() and
           constant() == 0;
}

void RoAffineExpression::set_multi_uncertainty_behaviour(
        RoAffineExpression::MultiUncertaintyBehaviour multi_uncertainty_behaviour) {
    _multi_uncertainty_behaviour = multi_uncertainty_behaviour;
}

RoAffineExpression::MultiUncertaintyBehaviour RoAffineExpression::multi_uncertainty_behaviour() const {
    return _multi_uncertainty_behaviour;
}

}