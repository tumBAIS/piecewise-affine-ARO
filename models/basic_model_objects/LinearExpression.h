//
// Created by simon on 14/06/2021.
//

#ifndef ROBUSTOPTIMIZATION_LINEAREXPRESSION_H
#define ROBUSTOPTIMIZATION_LINEAREXPRESSION_H

#include <string>
#include <vector>

#include "ScaledVariable.h"

namespace robust_model {

template<class V>
class LinearExpression {
public:
    LinearExpression() = default;

    explicit LinearExpression(std::vector<ScaledVariable<V>> const& scaled_variables);

    explicit LinearExpression(ScaledVariable<V> const& scaled_variable);

    explicit LinearExpression(V variable);

    LinearExpression<V>& operator*=(double scale);

    LinearExpression<V>& operator/=(double scale);

    LinearExpression<V> operator*(double scale) const;

    friend LinearExpression<V> operator*(double scale, LinearExpression<V> const& expr) {
        return expr * scale;
    }

    LinearExpression<V> operator/(double scale) const;

    LinearExpression<V>& operator+=(LinearExpression<V> const& other);

    LinearExpression<V>& operator+=(ScaledVariable<V> const& svar);

    LinearExpression<V>& operator+=(V const& var);

    template<LinearSubtractable<V> T>
    LinearExpression<V>& operator-=(T const& other);

    LinearExpression<V> operator-() const;

    template<LinearAddable<V> T>
    LinearExpression<V> operator+(T const& other) const;

    template<SolelyAffineAddable<V> T>
    AffineExpression<V> operator+(T const& other) const;

    friend AffineExpression<V> operator+(double other, LinearExpression<V> const& expr) {
        return AffineExpression<V>(other) + expr;
    }

    template<RobustVariableType V1, NonConstAffineConvertible<typename OtherRobustVariableType<V1>::Type::Reference> T>
    friend RoAffineExpression operator+(LinearExpression<VariableReference<V1>> const& v, T const& t);

    template<RobustVariableType V1, NonConstAffineConvertible<typename OtherRobustVariableType<V1>::Type::Reference> T>
    friend RoAffineExpression operator*(LinearExpression<VariableReference<V1>> const& v, T const& t);

    template<LinearSubtractable<V> T>
    LinearExpression<V> operator-(T const& other) const;

    template<SolelyAffineSubtractable<V> T>
    AffineExpression<V> operator-(T const& other) const;

    friend AffineExpression<V> operator-(double other, LinearExpression<V> const& expr) {
        return AffineExpression<V>(other) - expr;
    }

    template<RobustVariableType V1, NonConstAffineConvertible<typename OtherRobustVariableType<V1>::Type::Reference> T>
    friend RoAffineExpression operator-(LinearExpression<VariableReference<V1>> const& v, T const& t);

    template<AffineSubtractable<V> T>
    RawConstraint<AffineExpression<V>> operator>=(T const& other) const;

    template<AffineSubtractable<V> T>
    RawConstraint<AffineExpression<V>> operator<=(T const& other) const;

    template<AffineSubtractable<V> T>
    RawConstraint<AffineExpression<V>> operator==(T const& other) const;

    friend RawConstraint<AffineExpression<V>> operator>=(double other, LinearExpression<V> const& expr) {
        return other - expr >= 0;
    }

    friend RawConstraint<AffineExpression<V>> operator<=(double other, LinearExpression<V> const& expr) {
        return other - expr <= 0;
    }

    friend RawConstraint<AffineExpression<V>> operator==(double other, LinearExpression<V> const& expr) {
        return other - expr == 0;
    }

    template<RobustVariableType V1, NonConstAffineConvertible<typename OtherRobustVariableType<V1>::Type::Reference> T>
    friend RawConstraint<RoAffineExpression> operator>=(LinearExpression<VariableReference<V1>> const& v, T const& t);

    template<RobustVariableType V1, NonConstAffineConvertible<typename OtherRobustVariableType<V1>::Type::Reference> T>
    friend RawConstraint<RoAffineExpression> operator<=(LinearExpression<VariableReference<V1>> const& v, T const& t);

    template<RobustVariableType V1, NonConstAffineConvertible<typename OtherRobustVariableType<V1>::Type::Reference> T>
    friend RawConstraint<RoAffineExpression> operator==(LinearExpression<VariableReference<V1>> const& v, T const& t);

    template<class W>
    LinearExpression<W> translate_to_other(std::vector<W> const& new_vars) const;

    template<class W, AffineAddable<W> T>
    AffineExpression<W> substitute(std::vector<T> const& substitutions) const;

    std::vector<ScaledVariable<V>> const& scaled_variables() const;

    bool empty() const;

    std::string to_string() const;

    template<class S>
    double value(S const& solution) const;

    double lb() const;

    double ub() const;

public:
    template<class Iter>
    static LinearExpression<V> sum(Iter const& iter);

private:
    std::vector<ScaledVariable<V>> _scaled_variables;
};

}

#include "LinearExpression.tplt"

#endif //ROBUSTOPTIMIZATION_LINEAREXPRESSION_H
