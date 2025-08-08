import logging
from typing import List, Dict

import networkx as nx
from gseapy import enrichment_map
from pandas import DataFrame

from src.cellar3.analysis.functional.plots.plot import AbstractEnrichmentPlot
from src.cellar3.receipt import Receipt

logger = logging.getLogger('cellar.enrichment.plot.network')


class EnrichmentNetworkPlot(AbstractEnrichmentPlot):

    def __init__(self, top_terms=10, cutoff=0.05):
        self.cutoff = cutoff
        self.top_terms = top_terms

    def as_plot_data(self, df: DataFrame, meta: Dict[str, any], cutoff=0.05) -> Dict[str, List]:
        nodes, edges = enrichment_map(df, top_term=self.top_terms, cutoff=cutoff)
        logger.info(f"GSEapy Enrichment Map: {len(nodes)} nodes and {len(edges)} edges.")
        Receipt.add_gseapy_step(meta, name='extract nodes and edges', function='enrichment_map',
                                params={'top_term': self.top_terms, 'cutoff': cutoff})

        # Creates nodes connected by edges
        G = nx.from_pandas_edgelist(edges,
                                    source='src_idx',
                                    target='targ_idx',
                                    edge_attr=['jaccard_coef', 'overlap_coef', 'overlap_genes'])

        # Adds solitary nodes
        for idx in nodes.index:
            if idx not in G:
                G.add_node(idx)

        logger.info(f"NetworkX graph: {len(G.nodes)} nodes and {len(G.edges)} edges.")
        for idx, row in nodes.iterrows():
            G.nodes[idx].update(row.to_dict())

        positions = nx.spring_layout(G)
        edge_attributes = nx.get_edge_attributes(G, 'jaccard_coef')

        return {
            "nodes": [{
                'x': positions[node][0],
                'y': positions[node][1],
                'name': G.nodes[node]['Term'],
                'size': G.nodes[node]['Hits_ratio'],
                'color': G.nodes[node].get('NES') or G.nodes[node]['Hits_ratio']
            } for (node, row) in G.nodes(data=True)],
            "edges": [{
                'x': [positions[edge[0]][0], positions[edge[1]][0]],
                'y': [positions[edge[0]][1], positions[edge[1]][1]],
                'weight': edge_attributes.get(edge, 0)
            } for edge in G.edges()]
        }
