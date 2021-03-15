from itertools import chain

from collections import OrderedDict

from cobamp.core.models import ConstraintBasedModel
from cobamp.core.optimization import Solution
from cobamp.wrappers import get_model_reader
from cobamp.algorithms.kshortest import value_map_apply

import numpy as np

from imgrasp.graph.core import GenericGraph

class IntegratedGPRSolution(Solution):
    def  __init__(self, full_solution: Solution, reaction_mapping, gene_expression_prefix, **kwargs):
        ov = full_solution.objective_value()
        names, values = zip(*full_solution.var_values().items())
        metabolic_solution = OrderedDict()
        for k,v in reaction_mapping.items():
            metabolic_solution[names[k]] = values[v] if isinstance(v, int) else values[v[0]] - values[v[1]]

        gene_activity_solution = OrderedDict()
        for i,k in enumerate(names):
            if k.startswith(gene_expression_prefix):
                gene_activity_solution[k[len(gene_expression_prefix):]] = values[i]

        self.__gene_activity_solution = gene_activity_solution
        kwargs['names'] = list(metabolic_solution.keys())
        super().__init__(value_map={i:metabolic_solution[k] for i,k in enumerate(metabolic_solution.keys())},
                         status=full_solution.status(), objective_value=ov, **kwargs)


        self.remaining_reactions = {k:full_solution.var_values()[k] for k in set(names) -
                                    (set(metabolic_solution.keys()) | set(gene_activity_solution.keys()))}

    @property
    def gene_activity(self):
        return self.__gene_activity_solution


def get_linear_gpr_model(model_object, max_flux=1e4, gene_drain_prefix='GXD_', gene_metab_prefix='GXM_',
                         reaction_metab_prefix='RXM_', reaction_conv_prefix='RXD_', gpr_component_metab_prefix='GCOMP_',
                         and_component_prefix='COMP_', or_component_prefix='ORCOMP_'):

    model = get_model_reader(model_object)
    gpr_container = model.gene_protein_reaction_rules
    gene_names = gpr_container.get_genes()
    gene_names_map = dict(zip(gene_names, range(len(gene_names))))
    added_nodes = []
    rx_to_node = {}
    for r_id, i in zip(model.r_ids,range(len(gpr_container))):
        x = [frozenset([k for k in isoform if k in gene_names])
             for isoform in gpr_container.get_gpr_as_lists(i)]
        if len(x) > 0 and len(x[0]) > 0:
            for node_set in x:
                if node_set not in added_nodes:
                    added_nodes.append(node_set)
            rx_to_node[r_id] = x
    gene_node_dict = {i:frozenset([x]) for i, x in enumerate(gene_names)}
    intermediate_node_dict = {i:k for i, k in
                              enumerate(set(added_nodes) - set(gene_node_dict.values()))}

    node_union = {v:k for k,v in gene_node_dict.items()}
    node_union.update({v:k+len(gene_node_dict) for k,v in intermediate_node_dict.items()})
    reaction_dict = {r_id[0]:i  for i, r_id in enumerate(rx_to_node.items())}


    # G - gene by composite_node matrix (m genes by n extra_nodes)
    G = np.zeros([len(gene_node_dict), len(intermediate_node_dict)])
    int_nodes_list = list(intermediate_node_dict.values())
    for i,k in enumerate(int_nodes_list):
        g_idx = [gene_names_map[x] for x in k]
        G[g_idx,i] = -1

    # R - gene/intermediate on rows vs individual OR components for reactions
    # Q - identity: rows = reaction metabolites vs columns = OR components

    comp_rows, pseudo_identity = [], []
    row_dim = sum(map(len,[gene_node_dict,intermediate_node_dict]))
    for rx, nodes in rx_to_node.items():
        rows = np.zeros([row_dim,len(nodes)])
        for i,n in enumerate(nodes):
            rows[node_union[n],i] = -1
        comp_rows.append(rows)
        ident_row = np.zeros([len(reaction_dict),len(nodes)])
        ident_row[reaction_dict[rx],:] = 1
        pseudo_identity.append(ident_row)

    R, Q = (np.concatenate(x, axis=1) for x in [comp_rows, pseudo_identity])
    Rg, Rc = R[:len(gene_names),:], R[len(gene_names):,:]
    gn, cn, orcn, rxn = len(gene_names), len(intermediate_node_dict), R.shape[1], Q.shape[0]

    del R

    mat_blocks = [
        [np.eye(gn), G, Rg, np.zeros([gn, rxn])],
        [np.zeros([cn, gn]), np.eye(cn), Rc, np.zeros([cn, rxn])],
        [np.zeros([rxn, gn+cn]), Q, -np.eye(rxn)]
    ]

    rxd_ordered = {v:k for k,v in reaction_dict.items()}
    rxd_list_ordered = [rxd_ordered[i] for i in range(len(reaction_dict))]
    full_mat = np.vstack(list(map(np.hstack, mat_blocks)))

    del mat_blocks, Rg, Rc, R, Q

    rx_identifiers = [gene_drain_prefix+g for g in gene_names] + [and_component_prefix+str(i) for i in range(len(int_nodes_list))] + \
                     list(chain(*[[''.join([or_component_prefix,rx,'_'+str(i)]) for i in range(len(rx_to_node[rx]))]
                                  for rx in rxd_list_ordered])) + \
                     [reaction_conv_prefix+rx for rx in rxd_list_ordered]
    met_identifiers = [gene_metab_prefix+g for g in gene_names] + [gpr_component_metab_prefix+str(i) for i in range(len(int_nodes_list))] + \
                      [reaction_metab_prefix+r for r in rxd_list_ordered]

    gene_model = ConstraintBasedModel(S=full_mat,
                                      thermodynamic_constraints = [[0,max_flux] for i in range(full_mat.shape[1])],
                                      reaction_names=rx_identifiers, metabolite_names=met_identifiers)

    return gene_model, full_mat



def get_integrated_gpr_model(model: ConstraintBasedModel, max_flux=1e4, reaction_metab_prefix='RXM_',
                             reaction_conv_prefix='RXD_'):

    integrated_irrev_model, cb_model_rev_map = model.make_irreversible()
    M,N = [len(x) for x in [integrated_irrev_model.metabolite_names, integrated_irrev_model.reaction_names]]
    gpr_model, _ = get_linear_gpr_model(model, max_flux=max_flux)
    rx_to_add = [n for n in gpr_model.reaction_names if n.startswith(reaction_conv_prefix)]
    mt_to_add = [n for n in gpr_model.metabolite_names if n.startswith(reaction_metab_prefix)]

    mat_ident = np.zeros([len(rx_to_add), len(integrated_irrev_model.reaction_names)])
    for i, r in enumerate(rx_to_add):
        mat_ident[i, cb_model_rev_map[model.map_labels['reaction'][r[4:]]]] = -1

    n_gpr_met, n_irrcb_rx = len(gpr_model.metabolite_names), len(integrated_irrev_model.reaction_names)
    integrated_irrev_model.add_metabolites(np.zeros([n_gpr_met, n_irrcb_rx]), gpr_model.metabolite_names)
    integrated_irrev_model.set_stoichiometric_matrix(mat_ident, rows=mt_to_add, update_only_nonzero=True)

    del mat_ident

    padding = np.zeros([M, len(gpr_model.reaction_names)])
    int_model_mat = np.vstack([padding, gpr_model.get_stoichiometric_matrix()])

    del padding

    integrated_irrev_model.add_reactions(int_model_mat, gpr_model.bounds, gpr_model.reaction_names)
    integrated_irrev_model.remove_reactions(rx_to_add)
    return integrated_irrev_model, cb_model_rev_map, mt_to_add


def merge_linear_with_causal_model(model: ConstraintBasedModel, graph: GenericGraph, gene_metab_prefix='GXM_',
                                   gene_drain_prefix='GXD_', gene_pool_prefix='GXP_', gene_poolrx_prefix='ExGXP_',
                                   interaction_prefix='SRI_', max_flux=1e4):

    genes_in_causal = {g for g in graph.nodes}
    genes_in_model = {g[len(gene_metab_prefix):] for g in model.metabolite_names if g.startswith(gene_metab_prefix)}
    genes_to_add = genes_in_causal - genes_in_model
    total_genes = genes_in_causal | genes_in_model

    # remove all existing drains
    model.remove_reactions([gene_drain_prefix+g for g in genes_in_model])

    # add missing gene metabolites
    model.add_metabolites(
        np.zeros([len(genes_to_add), len(model.reaction_names)]),
        [gene_metab_prefix + x for x in genes_to_add])

    # add new pool metabolites
    model.add_metabolites(
        np.zeros([len(total_genes), len(model.reaction_names)]),
        [gene_pool_prefix + x for x in total_genes])

    # add reaction producing the gene metabolite and consuming its associated gene supply
    gxd_names, gxd_args = zip(*((gene_drain_prefix+g, {gene_metab_prefix+g: 1, gene_pool_prefix+g:-1}) for g in total_genes))
    gxp_names, gxp_args = zip(*((gene_poolrx_prefix+g, {gene_pool_prefix+g:1}) for g in total_genes))

    model.add_reactions(args=gxd_args+gxp_args, bounds=[(0, max_flux)]*(len(gxd_args)+len(gxp_args)), names=gxd_names+gxp_names)

    # get all of the graph's edges and store the source, target and associated weight
    graph_edges = graph.edges
    edge_dict = {interaction_prefix + i: (edge.source.identifier, edge.incident.identifier, edge.weight)
                 for i, edge in graph_edges.items()}

    sri_items = []
    for edge_name, edge_props in edge_dict.items():
        src, tar, weight = edge_props
        is_positive = weight >= 0
        edge_arg = {gene_metab_prefix + src: -1}
        if is_positive:
            edge_arg[gene_pool_prefix + tar] = 1
        else:
            edge_arg[gene_pool_prefix + tar] = -1
        sri_items.append([interaction_prefix+edge_name, edge_arg])

    sri_names, sri_args = zip(*sri_items)
    model.add_reactions(args=sri_args, bounds=[(0, max_flux)]*len(sri_args), names=sri_names)

    return edge_dict

class IntegratedGPRModel(ConstraintBasedModel):
    def __init__(self, model: ConstraintBasedModel, solver=None, gpr_rx_bound=1e4, **kwargs):

        integrated_irrev_model, cb_model_rev_map, mt_to_add = get_integrated_gpr_model(model, gpr_rx_bound, **kwargs)

        self.metabolic_fluxes = tuple(integrated_irrev_model.reaction_names)
        self.metabolic_rev_map = cb_model_rev_map

        S = integrated_irrev_model.get_stoichiometric_matrix()
        bounds = integrated_irrev_model.bounds
        rnames = integrated_irrev_model.reaction_names
        mnames = integrated_irrev_model.metabolite_names

        del integrated_irrev_model, cb_model_rev_map, mt_to_add

        super().__init__(S=S, thermodynamic_constraints=bounds, reaction_names=rnames, metabolite_names=mnames,
                         optimizer=True, solver=solver, gprs=None)

        ## TODO: variables that map metabolic reactions into their GPR control reactions

    def optimize(self, coef_dict=None, minimize=False):
        sol = super(IntegratedGPRModel, self).optimize(coef_dict, minimize)
        ## TODO: make prefix a parameter in the future
        return IntegratedGPRSolution(sol, self.metabolic_rev_map, 'GXD_',
                                     names=[self.reaction_names[k] for k in self.metabolic_rev_map])


class IntegratedGPRCausalModel(IntegratedGPRModel):
    def __init__(self, model: ConstraintBasedModel, graph: GenericGraph, solver=None, gpr_rx_bound=1e4,
                 gene_metab_prefix='GXM_',gene_drain_prefix='GXD_', gene_pool_prefix='GXP_',
                 gene_poolrx_prefix='ExGXP_', interaction_prefix='SRI_'):

        super().__init__(model, solver, gpr_rx_bound)

        self.__edges = merge_linear_with_causal_model(self, graph, gene_metab_prefix, gene_drain_prefix,
                                                     gene_pool_prefix, gene_poolrx_prefix, interaction_prefix)

    @property
    def causal_interaction_edges(self):
        return self.__edges