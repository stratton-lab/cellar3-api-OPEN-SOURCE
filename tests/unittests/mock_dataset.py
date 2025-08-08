from typing import Dict

from anndata import AnnData

from src.cellar3.datasets.dataset import Dataset
from tests.unittests.tools import AnnDataBuilder


def get_anndata_sample1():
    """
    Gene1: Has only one expressed cell in target, and none in background
    Gene2 : Has one expressed cell in target and one in background
    Gene3: Has one expressed cell in target
    :return:
    """
    return AnnDataBuilder(genes=['Gene1', 'Gene2', 'Gene3'],
                          cells=['Cell1_T', 'Cell2', 'Cell3_T', 'Cell4_B', 'Cell5', 'Cell6']) \
        .set_expressions('Gene1', {'Cell1_T': 0, 'Cell2': 1, 'Cell3_T': 0.5}) \
        .set_expressions('Gene2', {"Cell1_T": 1, 'Cell2': 1, 'Cell3_T': 1, 'Cell4_B': 0.05, 'Cell5': 1, 'Cell6': 1}) \
        .set_expressions('Gene3', {'Cell1_T': 1, 'Cell2': 2}) \
        .build()


def get_anndata_sample2():
    """
    Gene1: Has only one expressed cell in target, and none in background
    Gene2 : Has one expressed cell in target and one in background
    Gene3: Has one expressed cell in target
    :return:
    """
    return AnnDataBuilder(genes=['Gene1', 'Gene2', 'Gene3'],
                          cells=['Cell1_T', 'Cell2', 'Cell3_T', 'Cell4_B', 'Cell5_B', 'Cell6_B']) \
        .set_expressions('Gene1', {'Cell1_T': 0, 'Cell2': 1, 'Cell3_T': 0.5}) \
        .set_expressions('Gene2', {"Cell1_T": 1, 'Cell2': 1, 'Cell3_T': 1, 'Cell4_B': 1, 'Cell5_B': 1, 'Cell6_B': 1}) \
        .set_expressions('Gene3', {'Cell1_T': 1, 'Cell2': 2, 'Cell3_T': 10, 'Cell5_B': 0, 'Cell6_B': 20}) \
        .set_cell_obs('DATsubtype',
                      {'Cell1_T': 'Neuron', 'Cell5_B': 'Macrophage', 'Cell2': 'Neuron', 'Cell3_T': 'Neuron'}) \
        .set_2d_embedding('X_mock_embedding', [[0.1, 0.2],
                                               [0.3, 0.4],
                                               [0.5, 0.6],
                                               [0.7, 0.8],
                                               [0.9, 1.1],
                                               [1.2, 1.3]]) \
        .build()


def get_anndata_sample3():
    return AnnDataBuilder(genes=['BRAF', 'Igha', 'Igkc'],
                          cells=[f'Cell{i}_T' for i in range(1, 51)] + [f'Cell{i}_B' for i in range(1, 51)]) \
        .set_expressions('BRAF',
                         {**{f'Cell{i}_T': 1 for i in range(1, 51)}, **{f'Cell{i}_B': 0.1 for i in range(1, 51)}}) \
        .set_expressions('Igha',
                         {**{f'Cell{i}_T': 1 for i in range(1, 51)}, **{f'Cell{i}_B': 0.1 for i in range(1, 51)}}) \
        .set_expressions('Igkc',
                         {**{f'Cell{i}_T': 1 for i in range(1, 51)}, **{f'Cell{i}_B': 0.1 for i in range(1, 51)}}) \
        .build()


MOCK_DATASETS = {
    'example1': get_anndata_sample1(),
    'example2': get_anndata_sample2(),
    'example3': get_anndata_sample3()
}


class MockDataset(Dataset):

    def __init__(self, dataset_id: str, meta: Dict = None):
        super().__init__(dataset_id)
        if meta:
            self.meta = meta

    def _load_meta(self, dataset_id: str) -> Dict:
        return {
            'public': True
        }

    def _load(self, dataset_id: str) -> AnnData:
        return MOCK_DATASETS[dataset_id]
