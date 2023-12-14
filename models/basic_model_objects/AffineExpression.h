//
// Created by simon on 14/06/2021.
//

#ifndef ROBUSTOPTIMIZATION_AFFINEEXPRESSION_H
#define ROBUSTOPTIMIZATION_AFFINEEXPRESSION_H

#include <string>
#include <vector>

#include "LinearExpression.h"

namespace robust_model {

template<class V>
class AffineExpression {
public:
    AffineExpression() = default;

    AffineExpression(double constant, LinearExpression<V> const& linear);

    explicit AffineExpression(double constant);

    template<ConvertibleTo<LinearExpression<V> > T>
    explicit AffineExpression(T const& expr);

    AffineExpression<V>& operator*=(double scale);

    AffineExpression<V>& operator/=(double scale);

    AffineExpression<V> operator*(double scale) const;

    friend AffineExpression<V> operator*(double scale, AffineExpression<V> const& expr){
        return expr*scale;
    }

    AffineExpression<V> operator/(double scale) const;

    AffineExpression<V>& operator+=(AffineExpression<V> const& other);

    AffineExpression<V>& operator+=(LinearExpression<V> const& linear);

    AffineExpression<V>& operator+=(ScaledVariable<V> const& svar);

    AffineExpression<V>& operator+=(V const& var);

    AffineExpression<V>& operator+=(double constant);

    template<AffineSubtractable<V> T>
    AffineExpression<V>& operator-=(T const& other);

    AffineExpression<V> operator-() const;

    template<AffineAddable<V> T>
    AffineExpression<V> operator+(T const& other) const;

    friend AffineExpression<V> operator+(double other, AffineExpression<V> const& expr) {
        return AffineExpression<V>(other) + expr;
    }

    template<RobustVariableType V1, NonConstAffineConvertible<typename OtherRobustVariableType<V1>::Type::Reference> T>
    friend RoAffineExpression operator+(AffineExpression<VariableReference<V1>> const& v, T const& t);

    template<RobustVariableType V1, NonConstAffineConvertible<typename OtherRobustVariableType<V1>::Type::Reference> T>
    friend RoAffineExpression operator*(AffineExpression<VariableReference<V1>> const& v, T const& t);

    template<AffineSubtractable<V> T>
    AffineExpression<V> operator-(T const& other) const;

    friend AffineExpression<V> operator-(double other, AffineExpression<V> const& expr) {
        return AffineExpression<V>(other) - expr;
    }

    template<RobustVariableType V1, NonConstAffineConvertible<typename OtherRobustVariableType<V1>::Type::Reference> T>
    friend RoAffineExpression operator-(AffineExpression<VariableReference<V1>> const& v, T const& t);

    template<AffineSubtractable<V> T>
    RawConstraint<AffineExpression<V>> operator>=(T const& other) const;

    template<AffineSubtractable<V> T>
    RawConstraint<AffineExpression<V>> operator<=(T const& other) const;

    template<AffineSubtractable<V> T>
    RawConstraint<AffineExpression<V>> operator==(T const& other) const;

    friend RawConstraint<AffineExpression<V>> operator>=(double other, AffineExpression<V> const& expr) {
        return other - expr >= 0;
    }

    friend RawConstraint<AffineExpression<V>> operator<=(double other, AffineExpression<V> const& expr) {
        return other - expr <= 0;
    }

    friend RawConstraint<AffineExpression<V>> operator==(double other, AffineExpression<V> const& expr) {
        return other - expr == 0;
    }

    template<RobustVariableType V1, NonConstAffineConvertible<typename OtherRobustVariableType<V1>::Type::Reference> T>
    friend RawConstraint<RoAffineExpression> operator>=(AffineExpression<VariableReference<V1>> const& v, T const& t);

    template<RobustVariableType V1, NonConstAffineConvertible<typename OtherRobustVariableType<V1>::Type::Reference> T>
    friend RawConstraint<RoAffineExpression> operator<=(AffineExpression<VariableReference<V1>> const& v, T const& t);

    template<RobustVariableType V1, NonConstAffineConvertible<typename OtherRobustVariableType<V1>::Type::Reference> T>
    friend RawConstraint<RoAffineExpression> operator==(AffineExpression<VariableReference<V1>> const& v, T const& t);

    template<class W>
    AffineExpression<W> translate_to_other(std::vector<W> const& new_vars) const;

    template<class W, AffineAddable<W> T>
    AffineExpression<W> substitute(std::vector<T> const& substitutions) const;

    double constant() const;

    LinearExpression<V> const& linear() const;

    std::string to_string() const;

    template<class S>
    double value(S const& solution) const;

    double lb() const;
    double ub() const;

private:
    double _constant = 0;
    LinearExpression<V> _linear;
};

}

#include "AffineExpression.tplt"

#endif //ROBUSTOPTIMIZATION_AFFINEEXPRESSION_H
