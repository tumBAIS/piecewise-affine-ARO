//
// Created by simon on 14/06/2021.
//

#ifndef ROBUSTOPTIMIZATION_SOCEXPRESSION_H
#define ROBUSTOPTIMIZATION_SOCEXPRESSION_H

#include <vector>
#include <optional>

#include "../../helpers/helpers.h"
#include "types_and_constants.h"
#include "predeclarations_and_concepts.h"
#include "AffineExpression.h"
#include "NormedAffineVector.h"

namespace robust_model {

template<class V>
class SOCExpression;

template<class T, class V>
concept SOCAddable = requires(SOCExpression<V> e, T t) {
    e += t;
};

template<class T, class V>
concept SOCSubtractable = requires(SOCExpression<V> e, T t) {
    e += -t;
};

template<class V>
class SOCExpression {
    /* Implements an expression of the form ||Ax+b|| + cTx + d
     * The vector Ax+b is stored in _normed_vector
     * The affine part cTx + d is stored in _normed_vector
     */
public:
    SOCExpression() = default;

    template<ConvertibleTo<AffineExpression<typename V::Reference>> E>
    SOCExpression(NormedAffineVector<V> const& normed_vector,
                  E const& affine
    );

    template<ConvertibleTo<AffineExpression<typename V::Reference>> E>
    explicit SOCExpression(E const& affine);

    SOCExpression<V>& operator+=(SOCExpression<V> const& other);

    template<AffineAddable<typename V::Reference> T>
    SOCExpression<V>& operator+=(T const& other);

    template<SOCSubtractable<V> T>
    SOCExpression<V>& operator-=(T const& other);

    SOCExpression<V> operator-() const;

    template<SOCAddable<V> T>
    SOCExpression<V> operator+(T const& other) const;

    template<ConvertibleTo<AffineExpression<typename V::Reference>> T>
    friend SOCExpression<V> operator+(T const& left, SOCExpression<V> const& right) {
        return SOCExpression<V>(left) += right;
    }

    template<SOCSubtractable<V> T>
    SOCExpression<V> operator-(T const& other) const;

    template<ConvertibleTo<AffineExpression<typename V::Reference>> T>
    friend SOCExpression<V> operator-(T const& left, SOCExpression<V> const& right) {
        return SOCExpression<V>(left) -= right;
    }

    template<SOCSubtractable<V> T>
    RawConstraint<SOCExpression<V>> operator>=(T const& other) const;

    template<SOCSubtractable<V> T>
    RawConstraint<SOCExpression<V>> operator<=(T const& other) const;

    template<SOCSubtractable<V> T>
    RawConstraint<SOCExpression<V>> operator==(T const& other) const;

    template<ConvertibleTo<AffineExpression<typename V::Reference>> T>
    friend RawConstraint<SOCExpression<V>> operator>=(T const& other, RoAffineExpression const& expr) {
        return other - expr >= 0;
    }

    template<ConvertibleTo<AffineExpression<typename V::Reference>> T>
    friend RawConstraint<SOCExpression<V>> operator<=(T const& other, RoAffineExpression const& expr) {
        return other - expr <= 0;
    }

    template<ConvertibleTo<AffineExpression<typename V::Reference>> T>
    friend RawConstraint<SOCExpression<V>> operator==(T const& other, RoAffineExpression const& expr) {
        return other - expr == 0;
    }

    template<class W>
    SOCExpression<W> translate_to_other(std::vector<VariableReference<W>> const& new_vars) const;

    template<class W, AffineAddable<typename W::Reference> T>
    SOCExpression<W> substitute(std::vector<T> const& substitutions) const;

    bool is_affine() const;

    bool is_single_unscaled_variable() const;

    NormedAffineVector<V> const& normed_vector() const;

    AffineExpression<typename V::Reference> const& affine() const;

    std::string to_string() const;

    template<class S>
    double value(S const& solution) const;

    double lb() const;
    double ub() const;

public:
    template<class Iter>
    static SOCExpression<V> norm(Iter const& iter, VectorNormType norm_type=VectorNormType::Two);

private:
    std::optional<NormedAffineVector<V>> _normed_vector;
    AffineExpression<typename V::Reference> _affine;
};

}

#include "SOCExpression.tplt"

#endif //ROBUSTOPTIMIZATION_SOCEXPRESSION_H
