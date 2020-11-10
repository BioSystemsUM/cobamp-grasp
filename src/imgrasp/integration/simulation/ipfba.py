from typing import Dict, Iterable

from imgrasp.integration.models import IntegratedGPRCausalModel

ipfba_default_objective_weights = ((0.9,1.0),(1.0,1.5),(1.0,1.5))

class SimulationObjective(object):
    def __init__(self, coefficients: Dict[str, float], minimize: bool, sink_name: str):
        self.coefs = coefficients
        self.minimize = minimize
        self.sink_name = sink_name if len(coefficients) > 1 else list(coefficients.keys())[0]
        self.no_sink = len(self.coefs) == 1

    def get_bounds(self, model, lower_bound, upper_bound):
        if not self.no_sink:
            model.set_reaction_bounds('EX_'+self.sink_name, lb=0, ub=1e20)

        obj_sol = model.optimize(self.coefs, self.minimize)
        if obj_sol.status() == 'optimal':
            obj_val = obj_sol.objective_value()
            return (obj_val*x for x in [lower_bound, upper_bound])
        else:
            raise Exception('Model is infeasible')


def ipfba_default_objective_factory(growth_reaction, metabolic_reactions, non_consumed_sinks):
    growth = SimulationObjective({growth_reaction: 1}, False, '')
    parsimonious_enz = SimulationObjective({k: 1 for k in metabolic_reactions}, True, 'ipFBA_parsimonious_enz')
    parsimonious_sub = SimulationObjective({k: 1 for k in non_consumed_sinks}, True, 'ipFBA_parsimonious_sub')
    return [growth, parsimonious_enz, parsimonious_sub]


def constrain_model_by_prefix(model, d, prefix):
    for k,v in d.items():
        model.set_reaction_bounds(prefix+k, lb=v[0], ub=v[1])

def prepare_ipFBA_model(integrated_gpr_model, objective_configurations=(),
                          constrain_genes=None, constrain_expression=None,
                          gene_demand_prefix='ExGXP_', gene_expression_prefix='GXD_'):

    for d,p in zip([constrain_genes, constrain_expression],[gene_demand_prefix, gene_expression_prefix]):
        if d is not None: constrain_model_by_prefix(integrated_gpr_model, {k:(0, v) for k,v in d.items()}, p)

    for objclass, objperc in zip(objective_configurations, [[0, 1e4]]*len(objective_configurations)):
        if not objclass.no_sink:
            integrated_gpr_model.add_objective_constraint(objclass.coefs, objclass.minimize,
                                                          objperc, objclass.sink_name,
                                                          constrain=False)



def ipFBA(integrated_gpr_model, objectives: Iterable[SimulationObjective], objective_bounds, optimization_weights,
          expression_constraints=None, gene_constraints=None,
          gene_demand_prefix='ExGXP_', gene_expression_prefix='GXD_'):

    with integrated_gpr_model as m:
        for d,p in zip([gene_constraints, expression_constraints],[gene_demand_prefix, gene_expression_prefix]):
            if d is not None: constrain_model_by_prefix(integrated_gpr_model, {k:(0, v) for k,v in d.items()}, p)

        for objective, bounds in zip(objectives, objective_bounds):
            lb, ub = objective.get_bounds(m, *bounds)
            m.set_reaction_bounds('EX_'+objective.sink_name if not objective.no_sink
                                                     else list(objective.coefs.keys())[0], lb=lb, ub=ub)

        return m.optimize(optimization_weights, minimize=True)

class IPFBA(object):
    def __init__(self, icmodel: IntegratedGPRCausalModel, objectives: Iterable[SimulationObjective],
                 gene_demand_prefix='ExGXP_', gene_metab_prefix='GXM_', gene_expression_prefix='GXD_'):
        self.prefix = {'gdp': gene_demand_prefix, 'gmp':gene_metab_prefix, 'gep': gene_expression_prefix}
        self.model = icmodel
        self.objectives = objectives
        self.__remove_old_objective_constraints(objectives)

        prepare_ipFBA_model(integrated_gpr_model=self.model, objective_configurations=objectives,
                            gene_demand_prefix=gene_demand_prefix)

    def __remove_old_objective_constraints(self, objectives):
        to_remove = [ob.sink_name for ob in objectives if not ob.no_sink]
        reactions_to_remove = [k for k in ['EX_'+s for s in to_remove] if k in self.model.reaction_names]
        metabolites_to_remove = [k for k in to_remove if k in self.model.metabolite_names]
        if len(reactions_to_remove) > 0:
            print('Removing reactions:',reactions_to_remove)
            self.model.remove_reactions(reactions_to_remove)

        if len(metabolites_to_remove) > 0:
            print('Removing metabolites:',metabolites_to_remove)
            self.model.remove_metabolites(metabolites_to_remove)

    def run(self, objective_bounds, gene_constraints, expression_constraints, optimization_weights):
        return ipFBA(self.model, self.objectives, objective_bounds, expression_constraints=expression_constraints,
          gene_constraints=gene_constraints, optimization_weights=optimization_weights,
          gene_demand_prefix=self.prefix['gdp'], gene_expression_prefix=self.prefix['gep'])