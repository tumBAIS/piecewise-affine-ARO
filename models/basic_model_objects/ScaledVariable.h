//
// Created by simon on 14/06/2021.
//

#ifndef ROBUSTOPTIMIZATION_SCALEDVARIABLE_H
#define ROBUSTOPTIMIZATION_SCALEDVARIABLE_H

#include <string>

#include "VariableReference.h"

namespace robust_model {

template<class V>
class ScaledVariable {
public:
    ScaledVariable(double scale, V const& variable);

    explicit ScaledVariable(V const& variable);

    ScaledVariable<V>& operator*=(double scale);

    ScaledVariable<V>& operator/=(double scale);

    ScaledVariable<V> operator*(double scale) const;

    friend ScaledVariable<V> operator*(double scale, ScaledVariable<V> const& svar){
        return svar*scale;
    }

    ScaledVariable<V> operator/(double scale) const;

    ScaledVariable<V> operator-() const;

    template<LinearAddable<V> T>
    LinearExpression<V> operator+(T const& other) const;

    template<SolelyAffineAddable<V> T>
    AffineExpression<V> operator+(T const& other) const;

    friend AffineExpression<V> operator+(double other, ScaledVariable<V> v) {
        return AffineExpression<V>(other) + v;
    }

    template<RobustVariableType V1, NonConstAffineConvertible<typename OtherRobustVariableType<V1>::Type::Reference> T>
    friend RoAffineExpression operator+(ScaledVariable<VariableReference<V1> > const& v, T const& t);

    template<RobustVariableType V1, NonConstAffineConvertible<typename OtherRobustVariableType<V1>::Type::Reference> T>
    friend RoAffineExpression operator*(ScaledVariable<VariableReference<V1> > const& v, T const& t);

    template<LinearSubtractable<V> T>
    LinearExpression<V> operator-(T const& other) const;

    template<SolelyAffineSubtractable<V> T>
    AffineExpression<V> operator-(T const& other) const;

    friend AffineExpression<V> operator-(double other, ScaledVariable<V> v) {
        return AffineExpression<V>(other) - v;
    }

    template<RobustVariableType V1, NonConstAffineConvertible<typename OtherRobustVariableType<V1>::Type> T>
    friend RoAffineExpression operator-(ScaledVariable<typename V1::Reference> const& v, T const& t);

    template<AffineSubtractable<V> T>
    RawConstraint<AffineExpression<V>> operator>=(T const& other) const;

    template<AffineSubtractable<V> T>
    RawConstraint<AffineExpression<V>> operator<=(T const& other) const;

    template<AffineSubtractable<V> T>
    RawConstraint<AffineExpression<V>> operator==(T const& other) const;

    friend RawConstraint<AffineExpression<V>> operator>=(double other, ScaledVariable<V> const& v) {
        return other - v >= 0;
    }

    friend RawConstraint<AffineExpression<V>> operator<=(double other, ScaledVariable<V> const& v) {
        return other - v <= 0;
    }

    friend RawConstraint<AffineExpression<V>> operator==(double other, ScaledVariable<V> const& v) {
        return other - v == 0;
    }

    template<RobustVariableType V1, NonConstAffineConvertible<typename OtherRobustVariableType<V1>::Type> T>
    friend RawConstraint<RoAffineExpression> operator>=(ScaledVariable<typename V1::Reference> const& v, T const& t);

    template<RobustVariableType V1, NonConstAffineConvertible<typename OtherRobustVariableType<V1>::Type> T>
    friend RawConstraint<RoAffineExpression> operator<=(ScaledVariable<typename V1::Reference> const& v, T const& t);

    template<RobustVariableType V1, NonConstAffineConvertible<typename OtherRobustVariableType<V1>::Type> T>
    friend RawConstraint<RoAffineExpression> operator==(ScaledVariable<typename V1::Reference> const& v, T const& t);

    double scale() const;

    V const& variable() const;

    std::string to_string() const;

    template<class S>
    double value(S const& solution) const;

    double lb() const;
    double ub() const;

private:
    double _scale;
    V _variable;
};

}

#include "ScaledVariable.tplt"

#endif //ROBUSTOPTIMIZATION_SCALEDVARIABLE_H
