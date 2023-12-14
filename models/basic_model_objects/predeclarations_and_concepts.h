//
// Created by simon on 14/06/2021.
//

#ifndef ROBUSTOPTIMIZATION_PREDECLARATIONS_AND_CONCEPTS_H
#define ROBUSTOPTIMIZATION_PREDECLARATIONS_AND_CONCEPTS_H

#include <concepts>

namespace robust_model {

//*
// -----------
// Predeclarations
// -----------
// *//

template<class V>
class VariableReference;

template<class V>
class ScaledVariable;

template<class V>
class LinearExpression;

template<class V>
class AffineExpression;

template<class E>
class RawConstraint;

class DecisionVariable;

class UncertaintyVariable;

class RoAffineExpression;

//*
// -----------
// CONCEPTS
// -----------
// *//

template<class V>
concept RobustVariableType = std::same_as<V, DecisionVariable> || std::same_as<V, UncertaintyVariable>;

template<RobustVariableType V>
struct OtherRobustVariableType{
    using Type = void;
};

template<>
struct OtherRobustVariableType<DecisionVariable>{
    using Type = UncertaintyVariable;
};

template<>
struct OtherRobustVariableType<UncertaintyVariable>{
    using Type = DecisionVariable;
};


template<class T>
concept Negatable = requires(T t){
    -t;
};

template<class T, class E>
concept AddableToExpression = requires(E e, T t){
    e += t;
};

template<class T, class E>
concept SubtractableFromExpression = requires(E e, T t){
    e += -t;
};

template<class T, class V>
concept LinearAddable = AddableToExpression<T, LinearExpression<V>>;

template<class T, class V>
concept AffineAddable = AddableToExpression<T, AffineExpression<V>>;

template<class T, class V>
concept SolelyAffineAddable = AffineAddable<T, V> && !LinearAddable<T, V>;

template<class T, class V>
concept LinearSubtractable = SubtractableFromExpression<T, LinearExpression<V>>;

template<class T, class V>
concept AffineSubtractable = SubtractableFromExpression<T, AffineExpression<V>>;

template<class T, class V>
concept SolelyAffineSubtractable = AffineSubtractable<T, V> && !LinearSubtractable<T, V>;

template<class T, class V>
concept ConvertibleTo = requires(T t){ V(t); };

template <class T, class V>
concept NonConstAffineConvertible = ConvertibleTo<T, AffineExpression<V>> && !std::convertible_to<T, double>;

template<class S, class V>
concept SolutionProvider = requires(S s, V v){
    s.value(v);
};


}

#endif //ROBUSTOPTIMIZATION_PREDECLARATIONS_AND_CONCEPTS_H
