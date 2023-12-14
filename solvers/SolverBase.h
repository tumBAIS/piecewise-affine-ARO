//
// Created by thomae on 12/8/2023.
//

#ifndef ROBUSTOPTIMIZATION_SOLVERBASE_H
#define ROBUSTOPTIMIZATION_SOLVERBASE_H

#include <optional>
#include "../helpers/helpers.h"

namespace solvers {

class SolverBase {
public:
    enum class Status{
        UNSOLVED,
        OPTIMAL,
        TIME_LIMIT,
        MEMORY_LIMIT
    };

public:

    void build();

    bool built() const;

    void solve();
    // runtime limit in seconds
    void set_runtime_limit(double limit);
    void set_runtime_limit(std::optional<double> const& limit);

    // memory limit in GB
    void set_memory_limit(double limit);
    void set_memory_limit(std::optional<double> const& limit);

    Status status() const;

    bool has_solution() const;

    double runtime() const;

    double objective_value() const;

protected:
    void set_status(Status status);
    void set_runtime(double runtime);
    void set_objective_value(double objective_value);
    bool has_runtime_limit() const;
    double runtime_limit() const;
    std::optional<double> const& optional_runtime_limit() const;
    bool has_memory_limit() const;
    double memory_limit() const;
    std::optional<double> const& optional_memory_limit() const;

    void set_parameters_from_other(SolverBase const& other);
    void set_parameters_to_other(SolverBase & other) const;
    void set_results_from_other(SolverBase const& other);

private:
    virtual void build_implementation(){};
    virtual void update_implementation(){};
    virtual void solve_implementation() = 0;

private:
    Status _status = Status::UNSOLVED;
    std::optional<double> _runtime_limit;
    std::optional<double> _memory_limit;
    std::optional<double> _runtime;
    std::optional<double> _objective_value;
    bool _built = false;
    bool _encode_infeasible_results = true;
};

}

#endif //ROBUSTOPTIMIZATION_SOLVERBASE_H
