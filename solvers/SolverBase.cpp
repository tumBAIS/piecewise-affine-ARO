//
// Created by thomae on 12/8/2023.
//

#include "SolverBase.h"

namespace solvers {

void SolverBase::build() {
    helpers::exception_check(not built(), "Only build solver once!");
    build_implementation();
    _built = true;
}

bool SolverBase::built() const {
    return _built;
}

void SolverBase::solve() {
    if (not built()) {
        build();
    }
    update_implementation();
    solve_implementation();
}

void SolverBase::set_runtime_limit(double limit) {
    _runtime_limit = limit;
}

void SolverBase::set_memory_limit(double limit) {
    _memory_limit = limit;
}

bool SolverBase::has_runtime_limit() const {
    return _runtime_limit.has_value();
}

double SolverBase::runtime_limit() const {
    return _runtime_limit.value();
}

bool SolverBase::has_memory_limit() const {
    return _memory_limit.has_value();
}

double SolverBase::memory_limit() const {
    return _memory_limit.value();
}

void SolverBase::set_runtime_limit(std::optional<double> const& limit) {
    _runtime_limit = limit;
}

void SolverBase::set_memory_limit(std::optional<double> const& limit) {
    _memory_limit = limit;
}

std::optional<double> const& SolverBase::optional_runtime_limit() const {
    return _runtime_limit;
}

std::optional<double> const& SolverBase::optional_memory_limit() const {
    return _memory_limit;
}

void SolverBase::set_status(SolverBase::Status status) {
    _status = status;
    if (status == Status::UNSOLVED) {
        _runtime = {};
        _objective_value = {};
    }
}

SolverBase::Status SolverBase::status() const {
    return _status;
}

double SolverBase::runtime() const {
    if (_encode_infeasible_results and not has_solution())
        return std::nan("0");
    return _runtime.value();
}

bool SolverBase::has_solution() const {
    return status() == Status::OPTIMAL;
}

void SolverBase::set_runtime(double runtime) {
    _runtime = runtime;
}

double SolverBase::objective_value() const {
    if (_encode_infeasible_results and not has_solution())
        return std::nan("0");
    return _objective_value.value();
}

void SolverBase::set_objective_value(double objective_value) {
    _objective_value = objective_value;
}

void SolverBase::set_parameters_from_other(SolverBase const& other) {
    set_runtime_limit(other.optional_runtime_limit());
    set_memory_limit(other.optional_memory_limit());
}

void SolverBase::set_parameters_to_other(SolverBase& other) const {
    other.set_parameters_from_other(*this);
}

void SolverBase::set_results_from_other(SolverBase const& other) {
    set_status(other.status());
    if (other.has_solution()) {
        set_runtime(other.runtime());
        set_objective_value(other.objective_value());
    }
}

}
