//
// Created by simon on 14/06/2021.
//

#ifndef ROBUSTOPTIMIZATION_VARIABLEREFERENCE_H
#define ROBUSTOPTIMIZATION_VARIABLEREFERENCE_H

#include "predeclarations_and_concepts.h"

namespace robust_model {

template<class V>
class VariableReference {
public:
    explicit VariableReference(typename V::Index id);

    VariableReference(V const& var);

    ScaledVariable<typename V::Reference> operator-() const;

    ScaledVariable<typename V::Reference> operator*(double s) const;

    ScaledVariable<typename V::Reference> operator/(double s) const;

    friend ScaledVariable<typename V::Reference> operator*(double s, VariableReference<V> const& v) {
        return v * s;
    }

    template<LinearAddable<typename V::Reference> T>
    LinearExpression<typename V::Reference> operator+(T const& other) const;

    template<SolelyAffineAddable<typename V::Reference> T>
    AffineExpression<typename V::Reference> operator+(T const& other) const;

    friend AffineExpression<typename V::Reference> operator+(double other, VariableReference<V> v) {
        return AffineExpression<typename V::Reference>(other) + v;
    }

    template<RobustVariableType V1, NonConstAffineConvertible<typename OtherRobustVariableType<V1>::Type::Reference> T>
    friend RoAffineExpression operator+(VariableReference<V1> const& v, T const& t);

    template<RobustVariableType V1, NonConstAffineConvertible<typename OtherRobustVariableType<V1>::Type::Reference> T>
    friend RoAffineExpression operator*(VariableReference<V1> const& v, T const& t);

    template<LinearSubtractable<typename V::Reference> T>
    LinearExpression<typename V::Reference> operator-(T const& other) const;

    template<SolelyAffineSubtractable<typename V::Reference> T>
    AffineExpression<typename V::Reference> operator-(T const& other) const;

    friend AffineExpression<typename V::Reference> operator-(double other, VariableReference<V> v) {
        return AffineExpression<typename V::Reference>(other) - v;
    }

    template<RobustVariableType V1, NonConstAffineConvertible<typename OtherRobustVariableType<V1>::Type::Reference> T>
    friend RoAffineExpression operator-(VariableReference<V1> const& v, T const& t);

    template<AffineSubtractable<typename V::Reference> T>
    RawConstraint<AffineExpression<typename V::Reference>> operator>=(T const& other) const;

    template<AffineSubtractable<typename V::Reference> T>
    RawConstraint<AffineExpression<typename V::Reference>> operator<=(T const& other) const;

    template<AffineSubtractable<typename V::Reference> T>
    RawConstraint<AffineExpression<typename V::Reference>> operator==(T const& other) const;

    friend RawConstraint<AffineExpression<typename V::Reference>> operator>=(double other, VariableReference<V> const& v) {
        return other - v >= 0;
    }

    friend RawConstraint<AffineExpression<typename V::Reference>> operator<=(double other, VariableReference<V> const& v) {
        return other - v <= 0;
    }

    friend RawConstraint<AffineExpression<typename V::Reference>> operator==(double other, VariableReference<V> const& v) {
        return other - v == 0;
    }

    template<RobustVariableType V1, NonConstAffineConvertible<typename OtherRobustVariableType<V1>::Type::Reference> T>
    friend RawConstraint<RoAffineExpression> operator>=(VariableReference<V1> const& v, T const& t);

    template<RobustVariableType V1, NonConstAffineConvertible<typename OtherRobustVariableType<V1>::Type::Reference> T>
    friend RawConstraint<RoAffineExpression> operator<=(VariableReference<V1> const& v, T const& t);

    template<RobustVariableType V1, NonConstAffineConvertible<typename OtherRobustVariableType<V1>::Type::Reference> T>
    friend RawConstraint<RoAffineExpression> operator==(VariableReference<V1> const& v, T const& t);

    size_t raw_id() const { return _id.raw_id();}
    V const* operator->() const{return _id.operator->();}
    operator typename V::Index() const {return _id;}

    std::string to_string() const;

    template<SolutionProvider<typename V::Index> S>
    double value(S const& solution) const;

    double lb() const;
    double ub() const;

private:
    typename V::Index _id;
};

}

#include "VariableReference.tplt"

#endif //ROBUSTOPTIMIZATION_VARIABLEREFERENCE_H
