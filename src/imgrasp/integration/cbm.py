from cobamp.core.models import ConstraintBasedModel
from cobamp.wrappers.external_wrappers import get_model_reader

from itertools import chain
import numpy as np

def get_linear_gpr_model(model_object, max_flux=1e4):

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
    mat_blocks = [
        [np.eye(gn), G, Rg, np.zeros([gn, rxn])],
        [np.zeros([cn, gn]), np.eye(cn), Rc, np.zeros([cn, rxn])],
        [np.zeros([rxn, gn+cn]), Q, -np.eye(rxn)]
    ]
    rxd_ordered = {v:k for k,v in reaction_dict.items()}
    rxd_list_ordered = [rxd_ordered[i] for i in range(len(reaction_dict))]
    full_mat = np.vstack(list(map(np.hstack, mat_blocks)))

    print(full_mat.shape)
    rx_identifiers = ['GXD_'+g for g in gene_names] + ['COMP_'+str(i) for i in range(len(int_nodes_list))] + \
                     list(chain(*[['_'.join(['ORCOMP',rx,str(i)]) for i in range(len(rx_to_node[rx]))]
                                  for rx in rxd_list_ordered])) + \
                     ['RXD_'+rx for rx in rxd_list_ordered]
    met_identifiers = ['GXM_'+g for g in gene_names] + ['GCOMP_'+str(i) for i in range(len(int_nodes_list))] + \
                      ['RXM_'+r for r in rxd_list_ordered]

    gene_model = ConstraintBasedModel(S=full_mat,
                                      thermodynamic_constraints = [[0,max_flux] for i in range(full_mat.shape[1])],
                                      reaction_names=rx_identifiers, metabolite_names=met_identifiers)
    return gene_model, full_mat

def get_integrated_gpr_model(model: ConstraintBasedModel):
    cb_model_irrev, cb_model_rev_map = cb_model.make_irreversible()
    gpr_model, _ = get_linear_gpr_model(model, max_flux=10000)
    rx_to_add = [n for n in gpr_model.reaction_names if 'RXD_' in n]
    mt_to_add = [n for n in gpr_model.metabolite_names if 'RXM_' in n]

    mat_ident = np.zeros([len(rx_to_add), len(cb_model_irrev.reaction_names)])
    for i, r in enumerate(rx_to_add):
        mat_ident[i, cb_model_rev_map[cb_model.map_labels['reaction'][r[4:]]]] = -1

    integrated_irrev_model, _ = cb_model.make_irreversible()
    n_gpr_met, n_irrcb_rx = len(gpr_model.metabolite_names), len(integrated_irrev_model.reaction_names)
    integrated_irrev_model.add_metabolites(np.zeros([n_gpr_met, n_irrcb_rx]), gpr_model.metabolite_names)
    integrated_irrev_model.set_stoichiometric_matrix(mat_ident, rows=mt_to_add, update_only_nonzero=True)

    padding = np.zeros([len(cb_model_irrev.metabolite_names), len(gpr_model.reaction_names)])
    int_model_mat = np.vstack([padding, gpr_model.get_stoichiometric_matrix()])

    integrated_irrev_model.add_reactions(int_model_mat, gpr_model.bounds, gpr_model.reaction_names)
    integrated_irrev_model.remove_reactions(rx_to_add)
    return integrated_irrev_model, cb_model_rev_map, mt_to_add


if __name__ == '__main__':
    from urllib.request import urlretrieve
    from cobra.io import read_sbml_model

    path, _ = urlretrieve('https://github.com/SysBioChalmers/Human-GEM/raw/master/modelFiles/xml/HumanGEM.xml')
    model = read_sbml_model(path)

    cb_model = get_model_reader(model).to_cobamp_cbm('CPLEX')

    integrated_gpr_model, reversible_mapping, gene_metabolites = get_integrated_gpr_model(cb_model)

    integrated_gpr_model.optimize({'biomass_human': 1})

    #
    # from cobamp.wrappers.method_wrappers import KShortestEFMEnumeratorWrapper
    # algo = KShortestEFMEnumeratorWrapper(model=integrated_gpr_model, consumed=['glc__D_e'], produced=['succ_e'],
    #                               non_consumed=[], algorithm_type='kse_iterative', stop_criteria=10,
    #                               big_m_value=1000,big_m=True)
    #
    # sols = list(algo.get_enumerator())
    # import pandas as pd
    # df = pd.DataFrame(sols)