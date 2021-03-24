from typing import Dict, Iterable
from collections import Counter

import warnings
from numpy import zeros, ones
from imgrasp.integration.models import IntegratedGPRCausalModel

ipfba_default_objective_weights = ((0.9, 1.0), (1.0, 1.5), (1.0, 1.5))


class SimulationObjective(object):
    def __init__(self, coefficients: Dict[str, float], minimize: bool, sink_name: str):
        self.coefs = coefficients
        self.minimize = minimize
        self.sink_name = sink_name if len(coefficients) > 1 else list(coefficients.keys())[0]
        self.no_sink = len(self.coefs) == 1

    def get_bounds(self, model, lower_bound, upper_bound):
        if not self.no_sink:
            model.set_reaction_bounds('EX_' + self.sink_name, lb=0, ub=1e20)

        obj_sol = model.optimize(self.coefs, self.minimize)
        if obj_sol.status() == 'optimal':
            obj_val = obj_sol.objective_value()
            return (obj_val * x for x in [lower_bound, upper_bound])
        else:
            raise Exception('Model is infeasible')


def ipfba_default_objective_factory(growth_reaction, metabolic_reactions, non_consumed_sinks):
    growth = SimulationObjective({growth_reaction: 1}, False, '')
    parsimonious_enz = SimulationObjective({k: 1 for k in metabolic_reactions}, True, 'ipFBA_parsimonious_enz')
    parsimonious_sub = SimulationObjective({k: 1 for k in non_consumed_sinks}, True, 'ipFBA_parsimonious_sub')
    return [growth, parsimonious_enz, parsimonious_sub]


def constrain_model_by_prefix(model, d, prefix):
    for k, v in d.items():
        model.set_reaction_bounds(prefix + k, lb=v[0], ub=v[1])


def prepare_ipFBA_model(integrated_gpr_model, objective_configurations=(),
                        constrain_genes=None, constrain_expression=None,
                        gene_demand_prefix='ExGXP_', gene_expression_prefix='GXD_'):
    for d, p in zip([constrain_genes, constrain_expression], [gene_demand_prefix, gene_expression_prefix]):
        if d is not None: constrain_model_by_prefix(integrated_gpr_model, {k: (0, v) for k, v in d.items()}, p)

    for objclass, objperc in zip(objective_configurations, [[0, 1e4]] * len(objective_configurations)):
        if not objclass.no_sink:
            integrated_gpr_model.add_objective_constraint(objclass.coefs, objclass.minimize,
                                                          objperc, objclass.sink_name,
                                                          constrain=False)


def set_ipFBA_bounds(integrated_gpr_model, objectives: Iterable[SimulationObjective], objective_bounds,
                     expression_constraints=None, gene_constraints=None,
                     gene_demand_prefix='ExGXP_', gene_expression_prefix='GXD_'):
    for d, p in zip([gene_constraints, expression_constraints], [gene_demand_prefix, gene_expression_prefix]):
        if d is not None: constrain_model_by_prefix(integrated_gpr_model, {k: (0, v) for k, v in d.items()}, p)

    for objective, bounds in zip(objectives, objective_bounds):
        lb, ub = objective.get_bounds(integrated_gpr_model, *bounds)
        integrated_gpr_model.set_reaction_bounds('EX_' + objective.sink_name if not objective.no_sink
                                                 else list(objective.coefs.keys())[0], lb=lb, ub=ub)


class IPFBA(object):
    def __init__(self, icmodel: IntegratedGPRCausalModel, objectives: Iterable[SimulationObjective],
                 gene_demand_prefix='ExGXP_', gene_metab_prefix='GXM_', gene_expression_prefix='GXD_',
                 regulatory_interaction_prefix='SRI_'):
        self.prefix = {'gdp': gene_demand_prefix, 'gmp': gene_metab_prefix, 'gep': gene_expression_prefix,
                       'sri': regulatory_interaction_prefix}
        self.model = icmodel
        self.objectives = objectives
        self.__remove_old_objective_constraints(objectives)

        prepare_ipFBA_model(integrated_gpr_model=self.model, objective_configurations=objectives,
                            gene_demand_prefix=gene_demand_prefix)

    def __remove_old_objective_constraints(self, objectives):
        to_remove = [ob.sink_name for ob in objectives if not ob.no_sink]
        reactions_to_remove = [k for k in ['EX_' + s for s in to_remove] if k in self.model.reaction_names]
        metabolites_to_remove = [k for k in to_remove if k in self.model.metabolite_names]
        if len(reactions_to_remove) > 0:
            print('Removing reactions:', reactions_to_remove)
            self.model.remove_reactions(reactions_to_remove)

        if len(metabolites_to_remove) > 0:
            print('Removing metabolites:', metabolites_to_remove)
            self.model.remove_metabolites(metabolites_to_remove)

    def run(self, objective_bounds, gene_constraints, expression_constraints, optimization_weights):

        with self.model as m:
            set_ipFBA_bounds(m, self.objectives, objective_bounds,
                             expression_constraints=expression_constraints, gene_constraints=gene_constraints,
                             gene_demand_prefix=self.prefix['gdp'], gene_expression_prefix=self.prefix['gep'])

            return m.optimize(optimization_weights, minimize=True)

class LooplessIPFBA(IPFBA):
    def __init__(self, icmodel: IntegratedGPRCausalModel, objectives: Iterable[SimulationObjective], eps: float):

        super().__init__(icmodel, objectives)
        avars, acons, mmc_number = self.__apply_loopless_constraints(1, eps)
        self.__mmc_number = mmc_number

    def set_mismatch_relaxation(self, mismatch_relaxation):
        if self.__mmc_number > 0:
            ec = self.model.model.model.constraints['ipfba_mmsum']
            ec.lb = 0
            ec.ub = int(mismatch_relaxation*self.__mmc_number)
        else:
            raise Exception('Could not set mismatch relaxation - not enough mismatching interactions.')

    def __apply_loopless_constraints(self, mm_relaxation, eps):
        causal_edges = self.model.causal_interaction_edges
        gene_to_sri = {}

        for k, v in causal_edges.items():
            var_name = self.prefix['sri'] + k
            s, t, sg = v

            gene_exists = self.prefix['gmp'] + t in self.model.metabolite_names
            if t not in gene_to_sri.keys():
                gene_to_sri[t] = [[],[]]
            if gene_exists:
                gene_to_sri[t][0 if sg > 0 else 1].append(self.model.decode_index(var_name, 'reaction'))

        var_map = {}
        for k, v in gene_to_sri.items():
            for int_type, var_inds in zip(['activation', 'inhibition'], v):
                if len(var_inds) > 0:
                    var_map[k+'_'+int_type] = var_inds

        genes_to_constrain = [k for k,v in Counter(['_'.join(k.split('_')[:-1]) for k in var_map.keys()]).items() if v == 2]

        if len(genes_to_constrain) > 1:

            var_map = {k:v for k,v in var_map.items() if '_'.join(k.split('_')[:-1]) in genes_to_constrain}

            vars_to_add = list(var_map.keys())
            index_var_map = {v:k for k,v in enumerate(vars_to_add)}

            # i hate myself
            offset = len(self.model.model.model.variables)
            coffset = len(self.model.model.model.constraints)

            # add binary variables to represent active interaction types
            act_var_list = self.model.model.add_variables_to_model(var_names=vars_to_add,
                                                               lb=[0]*len(vars_to_add),
                                                               ub=[1]*len(vars_to_add),
                                                               var_types='binary')

            S_new = zeros([len(vars_to_add),offset+len(vars_to_add)])
            b_lb = [eps]*S_new.shape[0]
            b_ub = [None]*S_new.shape[0]
            b_zero = [0]*S_new.shape[0]

            for k,v in var_map.items():
                S_new[index_var_map[k], v] = 1

            irows = list(range(S_new.shape[0]))
            ivars = list(range(offset, offset+len(vars_to_add)))
            icomp_pos = list([1]*len(vars_to_add))
            icomp_neg = list([0]*len(vars_to_add))



            # for k,v in
            # # add constraints to activate each binary variable if its associated interaction type is active
            # for i in vars_to_add:
            #     var_map[]
            icpnames = ['looplessipfba_pos_'+str(i) for i in range(len(b_lb))]
            indcpos = self.model.model.add_rows_to_model(S_new, b_lb, b_ub, only_nonzero=True,
                                               indicator_rows=list(zip(irows, ivars, icomp_pos)), names=icpnames)

            icnnames = ['looplessipfba_neg_'+str(i) for i in range(len(b_lb))]
            indcneg = self.model.model.add_rows_to_model(S_new, b_zero, b_zero, only_nonzero=True,
                                               indicator_rows=list(zip(irows, ivars, icomp_neg)),names=icnnames)


            gvars_to_add = [g+'_mismatch' for g in genes_to_constrain]
            offset_ivm = len(index_var_map)
            index_var_map.update({v:k+offset_ivm for k,v in dict(enumerate(gvars_to_add)).items()})


            mm_var_list = self.model.model.add_variables_to_model(var_names=gvars_to_add,
                                                               lb=[0]*len(gvars_to_add),
                                                               ub=[1]*len(gvars_to_add),
                                                               var_types='binary')

            var_subset = act_var_list + mm_var_list

            S_mm = zeros([len(genes_to_constrain), len(var_subset)])

            mirows, mivars = list(range(S_mm.shape[0])), list(range(offset_ivm, offset_ivm+len(var_subset)))
            micomp_pos = list([1]*len(vars_to_add))
            micomp_neg = list([0]*len(vars_to_add))

            for i,k in enumerate(genes_to_constrain):
                S_mm[i,[index_var_map[k+suff] for suff in ['_activation','_inhibition']]] = 1

            icmpnames = ['looplessipfba_mmpos_'+str(i) for i in range(S_mm.shape[0])]
            mmcpos = self.model.model.add_rows_to_model(S_mm, [2]*S_mm.shape[0], [2]*S_mm.shape[0], only_nonzero=True,
                                               indicator_rows=list(zip(mirows, mivars, micomp_pos)), vars=var_subset,
                                                        names=icmpnames)

            icmnnames = ['looplessipfba_mmneg_'+str(i) for i in range(S_mm.shape[0])]
            mmcneg = self.model.model.add_rows_to_model(S_mm, [None]*S_mm.shape[0], [1]*S_mm.shape[0], only_nonzero=True,
                                               indicator_rows=list(zip(mirows, mivars, micomp_neg)), vars=var_subset,
                                                        names=icmnnames)

            mismatch_sum_mat = ones([1, len(mm_var_list)])

            sum_constraint = self.model.model.add_rows_to_model(
                mismatch_sum_mat, [0], [mm_relaxation*len(genes_to_constrain)], vars=mm_var_list, names=['ipfba_mmsum'])

            self.__sum_constraint = sum_constraint[0]

            return [v.name for v in var_subset], icpnames+icnnames+icmpnames+icmnnames+['ipfba_mmsum'], \
                   len(genes_to_constrain)
        else:
            warnings.warn('Could not set mismatch relaxation constraint - no mismatches are possible')
            return [], [], 0