//
// Created by simon on 18/06/2021.
//

#ifndef ROBUSTOPTIMIZATION_NORMEDAFFINEEXPRESSIONS_H
#define ROBUSTOPTIMIZATION_NORMEDAFFINEEXPRESSIONS_H

#include "AffineExpression.h"
#include "types_and_constants.h"
#include <numeric>
#include <functional>

namespace robust_model {

template<class V>
class NormedAffineVector {
public:
    NormedAffineVector(VectorNormType norm_type, std::vector<AffineExpression<typename V::Reference>> const& normed_vector);

    VectorNormType norm_type() const;

    std::vector<AffineExpression<typename V::Reference>> const& normed_vector() const;

    std::string to_string() const;

    template<class W>
    NormedAffineVector<W> translate_to_other(std::vector<VariableReference<W>> const& new_vars) const;

    template<class W, AffineAddable<typename W::Reference> T>
    NormedAffineVector<W> substitute(std::vector<T> const& substitutions) const;

    template<NonConstAffineConvertible<typename V::Reference> E>
    static NormedAffineVector<V> norm(std::vector<E> const& normed_vector,
                                      VectorNormType norm_type = VectorNormType::Two);

    static NormedAffineVector<V> norm(std::vector<AffineExpression<typename V::Reference>> const& normed_vector,
                                      VectorNormType norm_type = VectorNormType::Two);

    template<class S>
    double value(S const& solution) const;

    double lb() const;
    double ub() const;

    template<class W>
    requires requires(W w) {
        { w->id()} -> std::same_as<typename V::Index const&>;
    }
    bool is_variable_vector_norm(std::vector<W> const& variables) const;


private:
    double accumulated_vector_values(std::vector<double> const& values) const;

private:
    VectorNormType const _norm_type;
    std::vector<AffineExpression<typename V::Reference>> const _normed_vector;

};

}

#include "NormedAffineVector.tplt"

#endif //ROBUSTOPTIMIZATION_NORMEDAFFINEEXPRESSIONS_H
