//
// Created by simon on 12.05.21.
//

#ifndef ROBUSTOPTIMIZATION_DOMINATINGPOLICYSOLVER_H
#define ROBUSTOPTIMIZATION_DOMINATINGPOLICYSOLVER_H

#include "../../models/ROModel.h"
#include "../../models/SOCModel.h"
#include "../../helpers/helpers.h"
#include "AROPolicySolverBase.h"
#include "../soc_solvers/GurobiSOCSolver.h"

namespace robust_model {

class DominatingPolicySolver;

class Vertex: public helpers::IndexedObject<Vertex>{
public:
    using Index = typename helpers::IndexedObject<Vertex>::Index;

public:
    Vertex(Index id, DominatingPolicySolver const& dominating_model, std::vector<double> vertex);

    Vertex(Index const& id,
           DominatingPolicySolver const& dominating_model,
           std::vector<double> const& vertex,
           std::vector<UncertaintyVariable::Index> const& dominated_uncertainties);

    double vertex_value(UncertaintyVariable::Index const& var) const;

    void set_vertex_value(UncertaintyVariable::Index const& var, double value);

    std::vector<double> const& vertex_vector() const;

    void shift_vertex(UncertaintyVariable::Index uvar, double shift);

    std::vector<UncertaintyVariable::Index> const& dominated_uncertainties() const;

    std::string to_string() const;

private:
    DominatingPolicySolver const& _dominating_model;
    std::vector<double> _vertex;
    std::vector<UncertaintyVariable::Index> const _dominated_uncertainties;
};

class DominatingPolicySolver :
        public solvers::AROPolicySolverBase,
        private helpers::IndexedObjectOwner<Vertex>{
public:
    enum class Mode{
        BOX,
        BENTAL,
        THOMAE,
        THOMAE_SHIFT
    };

    using VertexId = typename Vertex::Index;
    using VertexOwner = helpers::IndexedObjectOwner<Vertex>;

    friend Vertex;
public:
    explicit DominatingPolicySolver(ROModel const& model, Mode mode);

    void set_reoptimize(bool reoptimize);

    std::size_t num_vertices() const;

    DominatingSolution dominating_solution(DecisionVariable::Index dvar) const;

    SolutionRealization specific_solution(std::vector<double> const& uncertainty_realization) const final;

    SolutionRealization specific_vertex_solution(std::size_t vertex) const;

    void build_vertices();

    void shift_vertices(UncertaintyVariable::Index uvar, double shift);

    VertexId dominating_vertex(UncertaintyVariable::Index uvar) const;

    VertexId vertex_id(std::size_t id) const;

    VertexId base_vertex() const;

    double rho_for_dominating_vertex(UncertaintyVariable::Index uvar) const;

    double approximation_factor() const;

    void merge_vertices(std::vector<std::vector<VertexId>> const& vertices);

private:

    void build_implementation() final;

    void solve_implementation() final;

    std::pair<double, double> get_beta_gamma(UncertaintySet::SpecialSetType set_type);

    std::pair<double, double> get_beta_gamma_small_approach(UncertaintySet::SpecialSetType set_type);

    VertexId add_vertex(std::vector<double> const& vertex,
                        std::vector<UncertaintyVariable::Index> const& dominated_variables);

    Vertex const& vertex(VertexId id) const;

    Vertex & vertex(VertexId id);

    void build_variables();

    void build_fix_period_constraints();

    void build_constraints();

    void build_objective();

    std::string first_failed_assumption() const;

    SOCVariable::Reference dominating_variable(DecisionVariable::Index const& var, VertexId vertex) const;

    SOCVariable::Reference uvar_vertex_shift(UncertaintyVariable::Index uvar) const;

    double vertex_value(UncertaintyVariable::Index const& var, VertexId vertex) const;

    double vertex_solution_value(UncertaintyVariable::Index const& var, VertexId vertex) const;

    AffineExpression<SOCVariable::Reference> vertex_expression(UncertaintyVariable::Index const& var, VertexId vertex) const;

    AffineExpression<SOCVariable::Reference> ro_affine_to_affine(RoAffineExpression const& old, VertexId v) const;

    double vertex_scale(VertexId vertex, std::vector<double> const& uncertainty_realization) const;

    double approach_scale() const;

    Vertex merge_vertices(std::vector<VertexId> vertices) const;

    void override_vertices(std::vector<Vertex> const& new_vertices);

    bool all_uncertainties_dominated_exactly_once();

    SOCModel & soc_model();
    SOCModel const& soc_model() const;

    solvers::SOCSolverBase & soc_solver();

    Mode const& mode() const;

    bool allow_shift() const;

private:
    SOCModel _soc_model;
    std::unique_ptr<solvers::GurobiSOCSolver> _soc_solver;
    std::vector<std::vector<SOCVariable::Reference>> _dominating_variables;
    Mode const _mode;
    bool _reoptimize = false;
    std::vector<SOCVariable::Reference> _vert_shift;
    std::vector<VertexId> _dominating_vertex_for_uvar;
};

}

#endif //ROBUSTOPTIMIZATION_DOMINATINGPOLICYSOLVER_H
