//
// Created by simon on 08.05.21.
//

#ifndef ROBUSTOPTIMIZATION_ROAFFINEEXPRESSION_H
#define ROBUSTOPTIMIZATION_ROAFFINEEXPRESSION_H

#include <vector>

#include "UncertaintyVariable.h"
#include "DecisionVariable.h"
#include "UncertaintyScaledDecision.h"
#include "AffineExpression.h"

#include "RawConstraint.h"

namespace robust_model {

class SolutionRealization;

template<class T>
concept NonCopyRoConvertible =ConvertibleTo<T, RoAffineExpression> && !std::same_as<T, RoAffineExpression>;

template<class T>
concept RoAddable = requires(RoAffineExpression r, T t) {
    r += t;
};

template<class T>
concept RoSubtractable = requires(RoAffineExpression r, T t) {
    r += -t;
};

class RoAffineExpression {
public:
    enum class MultiUncertaintyBehaviour{
        UNION,
        AVERAGE
    };

public:
    RoAffineExpression() = default;

    RoAffineExpression(double constant, LinearExpression<DecisionVariable::Reference> const& decisions,
                       LinearExpression<UncertaintyVariable::Reference> const& uncertainties,
                       LinearExpression<UncertaintyScaledDecision> const& uncertainty_decisions,
                       MultiUncertaintyBehaviour multi_uncertainty_behaviour);

    RoAffineExpression(double constant, LinearExpression<DecisionVariable::Reference> const& decisions,
                       LinearExpression<UncertaintyVariable::Reference> const& uncertainties,
                       LinearExpression<UncertaintyScaledDecision> const& uncertainty_decisions);

    RoAffineExpression(double constant,
                       LinearExpression<DecisionVariable::Reference> const& decisions,
                       LinearExpression<UncertaintyVariable::Reference> const& uncertainties);


    RoAffineExpression(double constant,
                       std::vector<ScaledVariable<DecisionVariable::Reference>> const& decisions,
                       std::vector<ScaledVariable<UncertaintyVariable::Reference>> const& uncertainties);

    template<ConvertibleTo<LinearExpression<DecisionVariable::Reference>> T>
    explicit RoAffineExpression(T const& t);

    template<ConvertibleTo<LinearExpression<UncertaintyVariable::Reference>> T>
    explicit RoAffineExpression(T const& t);

    template<ConvertibleTo<LinearExpression<UncertaintyScaledDecision>> T>
    explicit RoAffineExpression(T const& t);

    explicit RoAffineExpression(AffineExpression<DecisionVariable::Reference> const& t);

    explicit RoAffineExpression(AffineExpression<UncertaintyVariable::Reference> const& t);

    explicit RoAffineExpression(AffineExpression<UncertaintyScaledDecision> const& t);

    explicit RoAffineExpression(double c);

    RoAffineExpression& operator*=(double scale);

    RoAffineExpression& operator/=(double scale);

    RoAffineExpression operator*(double scale) const;

    friend RoAffineExpression operator*(double scale, RoAffineExpression const& expr) {
        return expr * scale;
    }

    RoAffineExpression operator/(double scale) const;

    RoAffineExpression& operator+=(RoAffineExpression const& other);

    template<LinearAddable<UncertaintyVariable::Reference> T>
    RoAffineExpression& operator+=(T const& expr);

    template<LinearAddable<DecisionVariable::Reference> T>
    RoAffineExpression& operator+=(T const& expr);

    template<LinearAddable<UncertaintyScaledDecision> T>
    RoAffineExpression& operator+=(T const& expr);

    RoAffineExpression& operator+=(double c);

    template<RobustVariableType V>
    RoAffineExpression& operator+=(AffineExpression<VariableReference<V> > const& expr);

    RoAffineExpression& operator+=(AffineExpression<UncertaintyScaledDecision> const& expr);

    RoAffineExpression operator-() const;

    template<RoSubtractable T>
    RoAffineExpression& operator-=(T const& t);

    template<RoAddable T>
    RoAffineExpression operator+(T const& t) const;

    template<NonCopyRoConvertible T>
    friend RoAffineExpression operator+(T const& other, RoAffineExpression const& expr) {
        return RoAffineExpression(other) += expr;
    }

    template<RoSubtractable T>
    RoAffineExpression operator-(T const& t) const;

    template<NonCopyRoConvertible T>
    friend RoAffineExpression operator-(T const& other, RoAffineExpression const& expr) {
        return RoAffineExpression(other) -= expr;
    }

    template<RoSubtractable T>
    RawConstraint<RoAffineExpression> operator>=(T const& other) const;

    template<RoSubtractable T>
    RawConstraint<RoAffineExpression> operator<=(T const& other) const;

    template<RoSubtractable T>
    RawConstraint<RoAffineExpression> operator==(T const& other) const;

    template<NonCopyRoConvertible T>
    friend RawConstraint<RoAffineExpression> operator>=(T const& other, RoAffineExpression const& expr) {
        return other - expr >= 0;
    }

    template<NonCopyRoConvertible T>
    friend RawConstraint<RoAffineExpression> operator<=(T const& other, RoAffineExpression const& expr) {
        return other - expr <= 0;
    }

    template<NonCopyRoConvertible T>
    friend RawConstraint<RoAffineExpression> operator==(T const& other, RoAffineExpression const& expr) {
        return other - expr == 0;
    }

    template<AffineAddable<typename DecisionVariable::Reference> TD, AffineAddable<typename UncertaintyVariable::Reference> TU>
    RoAffineExpression substitute(std::vector<TD> const& decision_substitutions,
                                  std::vector<TU> const& uncertainty_substitutions) const;

    template<class W, AffineAddable<W> TD, AffineAddable<W> TU>
    AffineExpression<W> substitute_to_other_affine_no_uncertain_decisions(std::vector<TD> const& decision_substitutions,
                                                   std::vector<TU> const& uncertainty_substitutions) const;

    template<class W, AffineAddable<W> TD, AffineAddable<W> TU>
    AffineExpression<W> substitute_to_other_affine(std::vector<TD> const& decision_substitutions,
                                                   std::vector<TU> const& uncertainty_substitutions) const;

    double constant() const;

    LinearExpression<DecisionVariable::Reference> const& decisions() const;

    LinearExpression<UncertaintyVariable::Reference> const& uncertainties() const;

    LinearExpression<UncertaintyScaledDecision> const& uncertainty_decisions() const;

    std::string to_string() const;

    double value(SolutionRealization const& realization) const;

    double lb() const;

    double ub() const;

    bool is_zero() const;

    void set_multi_uncertainty_behaviour(MultiUncertaintyBehaviour multi_uncertainty_behaviour);

    MultiUncertaintyBehaviour multi_uncertainty_behaviour() const;

private:
    double _constant = 0;
    LinearExpression<DecisionVariable::Reference> _decisions;
    LinearExpression<UncertaintyVariable::Reference> _uncertainties;
    LinearExpression<UncertaintyScaledDecision> _uncertainty_decisions;
    MultiUncertaintyBehaviour _multi_uncertainty_behaviour = MultiUncertaintyBehaviour::UNION;
};

}

#include "RoAffineExpression.tplt"

#endif //ROBUSTOPTIMIZATION_ROAFFINEEXPRESSION_H
