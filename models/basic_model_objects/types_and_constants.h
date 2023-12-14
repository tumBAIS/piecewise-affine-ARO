//
// Created by simon on 08/06/2021.
//

#ifndef ROBUSTOPTIMIZATION_TYPES_AND_CONSTANTS_H
#define ROBUSTOPTIMIZATION_TYPES_AND_CONSTANTS_H

#include <string>
#include "../../helpers/helpers.h"
#include "gurobi_c++.h"

namespace robust_model{
static constexpr double NO_VARIABLE_UB = std::numeric_limits<double>::infinity();
static constexpr double NO_VARIABLE_LB = -std::numeric_limits<double>::infinity();

using period_id = int;

enum class ConstraintSense {
    GEQ, LEQ, EQ
};

static std::string sense_string(ConstraintSense s){
    switch (s) {
        case ConstraintSense::GEQ:
            return ">=";
        case ConstraintSense::LEQ:
            return "<=";
        case ConstraintSense::EQ:
            return "==";
    }
    return "";
}

static char to_grb_sense(ConstraintSense s){
    switch (s) {
        case ConstraintSense::GEQ:
            return GRB_GREATER_EQUAL;
        case ConstraintSense::LEQ:
            return GRB_LESS_EQUAL;
        case ConstraintSense::EQ:
            return GRB_EQUAL;
    }
    return char();
}

enum class ObjectiveSense {
    MIN, MAX
};

static std::string sense_string(ObjectiveSense sense) {
    switch (sense) {
        case ObjectiveSense::MIN:
            return "MIN";
        case ObjectiveSense::MAX:
            return "MAX";
        default:
            helpers::exception_check(false, "Forbidden Case!");
            return "";
    }
}

static int to_grb_sense(ObjectiveSense sense){
    switch (sense) {
        case ObjectiveSense::MIN:
            return GRB_MINIMIZE;
        case ObjectiveSense::MAX:
            return GRB_MAXIMIZE;
        default:
            helpers::exception_check(false, "Forbidden Case!");
            return 0;
    }
}

enum class VectorNormType{
    One, Two, Max
};

static VectorNormType dual_norm_type(VectorNormType norm_type){
    switch (norm_type) {
        case VectorNormType::One:
            return VectorNormType::Max;
        case VectorNormType::Two:
            return VectorNormType::Two;
        case VectorNormType::Max:
            return VectorNormType::One;
        default:
            helpers::exception_check(false, "Forbidden Case!");
            return VectorNormType::One;
    }
}

static std::string type_string(VectorNormType norm_type){
    switch (norm_type) {
        case VectorNormType::One:
            return "One";
        case VectorNormType::Two:
            return "Two";
        case VectorNormType::Max:
            return "Max";
        default:
            helpers::exception_check(false, "Forbidden Case!");
            return "";
    }
}

enum class VariableType{
    Continuous, Binary, Integer
};

static std::string type_string(VariableType var_type){
    switch (var_type) {
        case VariableType::Continuous:
            return "Continuous";
        case VariableType::Binary:
            return "Binary";
        case VariableType::Integer:
            return "Integer";
        default:
            helpers::exception_check(false, "Forbidden Case!");
            return "";
    }
}

static char to_grb_type(VariableType sense){
    switch (sense) {
        case VariableType::Continuous:
            return GRB_CONTINUOUS;
        case VariableType::Binary:
            return GRB_BINARY;
        case VariableType::Integer:
            return GRB_INTEGER;
        default:
            helpers::exception_check(false, "Forbidden Case!");
            return 0;
    }
}

}


#endif //ROBUSTOPTIMIZATION_TYPES_AND_CONSTANTS_H
