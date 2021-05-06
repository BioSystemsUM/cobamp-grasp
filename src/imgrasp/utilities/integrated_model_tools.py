import pandas as pd

from itertools import chain
from typing import Tuple, Dict

from imgrasp.integration.models import IntegratedGPRModel, IntegratedGPRCausalModel

GraphInfoDict = Dict[str, object]


def integrated_model_to_graph(integrated_model: IntegratedGPRModel) -> Tuple[GraphInfoDict, GraphInfoDict]:
    """
    This function converts an IntegratedGPRModel class into its graph representation.
    :param integrated_model: An IntegratedGPRModel instance (or subclass)
    :return: A tuple of dictionaries - The first element maps edge names with its edge properties while the second
    tuple does the same for nodes

    The fields included in each of the dictionaries' values are the following:
    - name: A unique identifier for the entity
    - source (edges only): The source node from which the edge begins
    - target (edges only): The target node where the edge ends
    - igr_reaction: Corresponding reaction in the integrated model,
    - value: empty field (should be populated with a solution flux value for the igr_reaction,
    - type: one of the following for edges: (consumed, expression, activation, inhibition, produces)
            one of the following for nodes: (gene, reaction, metabolite)
    """

    if isinstance(integrated_model, IntegratedGPRCausalModel):
        causal_ndict = {name: {'name': name,
                               'source': edge_properties[0], 'target': edge_properties[1],
                               'igr_reaction': 'SRI_' + name, 'value': None,
                               'type': 'activation' if edge_properties[2] > 0 else 'inhibition'}
                        for name, edge_properties in integrated_model.causal_interaction_edges.items()}
    else:
        causal_ndict = {}

    cbmat = integrated_model.get_stoichiometric_matrix(rows=integrated_model.metabolic_compound_names,
                                                       columns=integrated_model.metabolic_fluxes)
    import scipy.sparse as sprs

    sprs_cbm_mat = sprs.csr_matrix(pd.DataFrame(cbmat,
                                                columns=integrated_model.metabolic_fluxes,
                                                index=integrated_model.metabolic_compound_names)).todok()

    gpr_ndict = {}
    for rx, maps in integrated_model.gpr_mapping.items():
        for gene_sets, igr_rx in zip(*maps):
            for g in gene_sets:
                name = igr_rx + '_' + g
                gpr_ndict[name] = {'name': name, 'source': g, 'target': rx, 'value': None,
                                   'igr_reaction': igr_rx, 'type': 'expression'}

    metabolic_ndict = {'MTBL' + str(k).replace(' ', ''): {
        'name': str(k),
        'source': integrated_model.metabolic_fluxes[k[1]] if v > 0 else integrated_model.metabolic_compound_names[k[0]],
        'target': integrated_model.metabolic_compound_names[k[0]] if v > 0 else integrated_model.metabolic_fluxes[k[1]],
        'value': None,
        'igr_reaction': integrated_model.metabolic_fluxes[k[1]], 'type': 'produced' if v > 0 else 'consumed'}
        for k, v in sprs_cbm_mat.items()}

    full_dict = {}
    for d in [metabolic_ndict, causal_ndict, gpr_ndict]:
        full_dict.update(d)

    rx_node_dict = {k: {'name': k, 'igr_reaction': k, 'value': None, 'type': 'reaction'}
                    for k in integrated_model.metabolic_fluxes}

    met_node_dict = {k: {'name': k, 'igr_reaction': None, 'value': None, 'type': 'metabolite'}
                     for k in integrated_model.metabolic_compound_names}

    genes = {v['source'] for k, v in gpr_ndict.items()} | set(
        chain(*[[v['source'], v['target']] for v in causal_ndict.values()]))

    gene_node_dict = {k: {'name': k, 'igr_reaction': 'GXD_' + k, 'value': None, 'type': 'gene'} for k in genes}

    full_node_dict = {}
    for d in [rx_node_dict, met_node_dict, gene_node_dict]:
        full_node_dict.update(d)

    return full_dict, full_node_dict
