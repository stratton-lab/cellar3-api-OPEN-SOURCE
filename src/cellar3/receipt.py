import datetime
import sys
from typing import Dict, List

import gseapy
import scanpy


class Receipt:

    def __init__(self):
        pass

    @classmethod
    def create(cls) -> Dict:
        return {
            "date": str(datetime.datetime.now()),
            "environment": {"python": str(sys.version)},
            "steps": [],
            "preprocessing": [],
            "disclaimer": "Analysis produced by this website will not produce identical results to original "
                          "manuscript due to real-time constraints. To reproduce original data, use GEO and GitHub "
                          "resources from original publications."
        }

    @classmethod
    def add_step(cls, receipt: Dict, step: any):
        if receipt:
            receipt.setdefault('steps', []).append(step)

    @classmethod
    def add_dataset_load_step(cls, receipt: Dict, dataset_meta: Dict):
        if receipt:
            Receipt.add_step(receipt, {
                'step': 'load dataset',
                'datasetId': dataset_meta['id']
            })

    @classmethod
    def add_scanpy_step(cls, receipt: Dict, name: str, function: str, params: Dict):
        if receipt:
            Receipt.add_step(receipt, {
                'step': name,
                'package': {'name': 'scanpy', 'version': scanpy.__version__},
                'function': function,
                'params': params
            })

    @classmethod
    def add_scanpy_subsample_step(cls, receipt: Dict, n_obs: int):
        Receipt.add_scanpy_step(receipt, name='subsampling',
                                function='sc.pp.subsample', params={'n_obs': n_obs})

    @classmethod
    def add_custom_subsample_step(cls, receipt: Dict, group_name: str, group_max_cells: int,
                                  max_cells_per_cluster: int,
                                  before_clusters: Dict[str, List[str]],
                                  after_clusters: Dict[str, List[str]],
                                  before_cells: List[str], after_cells: List[str]):
        if receipt:
            Receipt.add_step(receipt, {
                'step': f'subsampling {group_name}',
                'description': 'Sampling the cell ids into uniform cluster representation by cell types.',
                'groupMaxCells': group_max_cells,
                'maxCellsPerCluster': max_cells_per_cluster,
                'beforeClusters': {k: len(v) for k, v in before_clusters.items()},
                'afterClusters': {k: len(v) for k, v in after_clusters.items()},
                'beforeCells': len(before_cells),
                'afterCells': len(after_cells)
            })

    @classmethod
    def add_hvg_step(cls, receipt: Dict):
        Receipt.add_scanpy_step(receipt, name='filtering by highly variable genes',
                                function='sc.pp.highly_variable_genes', params={})

    @classmethod
    def add_dge_step(cls, receipt: Dict, method: str):
        Receipt.add_scanpy_step(receipt, name='differential gene expression',
                                function='sc.tl.rank_genes_groups', params={'method': method})

    @classmethod
    def add_gseapy_step(cls, receipt: Dict, name: str, function: str, params: Dict):
        if receipt:
            Receipt.add_step(receipt, {
                'step': name,
                'package': {'name': 'gseapy', 'version': gseapy.__version__},
                'function': function,
                'params': params
            })

    @classmethod
    def add_enrichir_step(cls, receipt: Dict, params: Dict):
        Receipt.add_gseapy_step(receipt, name='functional enrichment', function='gp.enrichr', params=params)

    @classmethod
    def add_gene_filter_step(cls, receipt: Dict, min_cells_group: int):
        if receipt:
            Receipt.add_step(receipt, {
                'step': 'filter genes by expression',
                'description': 'Only keeps genes that are expressed in more that n cells in both target and background.',
                'min_cells_group': min_cells_group
            })

    @classmethod
    def assign_category(cls, receipt: Dict, params: Dict):
        if receipt:
            Receipt.add_step(receipt, {
                'step': 'assign category',
                'thresholds': params
            })

    @classmethod
    def bar_plot(cls, receipt: Dict, params: Dict):
        if receipt:
            Receipt.add_step(receipt, {
                'step': 'Bar Plot',
                'params': params
            })
