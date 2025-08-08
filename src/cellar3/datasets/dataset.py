import logging
import math
import os
import random
from functools import lru_cache
from typing import Dict, List, Tuple, Union

import anndata
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objs as go
import scanpy as sc
from anndata import AnnData
from anndata._core.views import ArrayView
from django.conf import settings
from numpy import ndarray
from pandas import Series

from src.cellar3.datasets.meta_reader import META_READER
from src.cellar3.exceptions import UserException
from src.cellar3.receipt import Receipt

logger = logging.getLogger('cellar.dataset')

random.seed(42)


class DatasetException(Exception):
    pass


class Dataset:
    """
    Based on: https://github.com/euxhenh/cellar/blob/main/controller/cellar/core/_plots.py#L49
    Notes on X :
    - It's a CSCDataset : Compressed Sparse Column: columns (genes) are stored efficiently
    - Accessing: X[cell_idx, gene_idx]
    """

    def __init__(self, dataset_id: str):
        """
        :param dataset_id:
        Returning both public and private datasets.
        """
        self.dataset_id = dataset_id
        self.meta = self._load_meta(dataset_id)
        self.data = self._load(dataset_id)

    def _load_meta(self, dataset_id: str) -> Dict:
        return META_READER.get_dataset_info(dataset_id)

    @staticmethod
    @lru_cache(maxsize=20)
    def __read_file(path: str) -> AnnData:
        return anndata.read_h5ad(path, 'r')

    def _load(self, dataset_id: str) -> AnnData:
        """
        Based on:
        1. https://github.com/euxhenh/cellar/blob/main/controller/data_loader.py
        2. https://github.com/euxhenh/cellar/blob/a1b9e551ed1ba29a7808df9dc9cc4359aae27008/controller/cellar/core/_tools.py#L10
        In AnnData :
            - obsm stores 2D embeddings used to draw the plot
            - obs stored the labels for cells
        :param dataset_id:
        :return:
        """
        if not self.meta:
            raise DatasetException(f'No meta data loaded for dataset {dataset_id}.')

        file_path = self.meta['file']
        if not file_path:
            raise DatasetException(f'No file path declared for dataset {dataset_id}')

        path = f'{settings.DATASETS_FOLDER}/{file_path}'
        if not os.path.exists(path):
            logger.warning(f'Dataset not found for file path: {path}')
            raise DatasetException(f'Dataset {dataset_id} not available.')

        data = self.__read_file(path)
        self._qc(data)
        logger.info(f'Dataset contains {data.var.index.size} genes and {data.obs.index.size} cells. '
                    f'2D Embeddings: {list(data.obsm.keys())}')
        # logger.info(f'Available group columns: {list(data.obs.keys())}')

        return data

    def get_dataset_info(self):
        return {
            'embeddings': list(self.data.obsm.keys()),
            'genes_count': self.data.var.index.size,
            'cells_count': self.data.obs.index.size,
            'group_columns': list(self.data.obs.keys())
        }

    def get_genes_count(self):
        return self.data.var.index.size

    def get_genes(self) -> List[str]:
        """
        Returns the list of gene names in the dataset.
        :return:
        """
        return self.data.var_names.tolist()

    def get_cells_count(self):
        return self.data.obs.index.size

    def get_interactions(self, source: str, targets: List[str], interaction_type: Union[str, None]) -> Dict:
        """
        Returns a list of interactions between a single source and multiple targets (Cell Types)
        Optionally, interactions can be filtered by type.
        :param source:
        :param targets:
        :param interaction_type:
        :return:
        """
        if not self.has_interactions():
            raise UserException('Interaction data not available for this dataset.')

        df = self.data.uns['interactions']
        filtered_df = df[(df['source'] == source) & (df['target'].isin(targets))]
        # Available types are not restricted to type filtering
        available_types = filtered_df['annotation'].unique().tolist()
        if interaction_type is not None:
            filtered_df = filtered_df[filtered_df['annotation'] == interaction_type]
        return {
            'interactions': filtered_df.to_dict(orient='records'),
            'types': available_types
        }

    def has_interactions(self) -> bool:
        """
        Returns true if the dataset contains interactions.
        :return:
        """
        return 'interactions' in self.data.uns

    def get_available_embeddings(self):
        return list(self.data.obsm.keys())

    def _qc(self, data):
        if not data.var.index.is_unique:
            raise DatasetException('Gene names are not unique.')
        if data.shape[0] <= 2:
            raise DatasetException('Dataset contains less than 3 points.')
        if not data.obsm.keys():
            raise DatasetException('No Embeddings found')

        if self.meta:
            declared_cells_count = self.meta.get('cells')
            if data.n_obs != declared_cells_count:
                raise UserException(f'Incorrect number of cells: actual {data.n_obs}, declared {declared_cells_count}')

    def as_plot(self, embedding_name: str, group_by: str) -> go.Figure:
        """

        :param group_by: Column with 2D coordinates of cells
        :param embedding_name: Will use default embedding if None
        @todo Use Categorical order in 'color'
        :return:
        """
        logger.info(f'Building figure from embedding {embedding_name}, grouped by {group_by}')
        x, y = self.get_embedding(embedding_name=embedding_name)
        return px.scatter(
            x=x,
            y=y,
            color=self.data.obs[group_by].astype('str') if group_by in self.data.obs else None,  # Clusters
            custom_data=self.get_tooltip_data()
        )

    def as_plot_data(self, embedding_name: str, group_by: str) -> List[Dict]:
        """
        Returns a JSON representation, with a list of dictionaries, one for each trace.
        :param embedding_name:
        :param group_by:
        :return:
        """
        plot = self.as_plot(embedding_name, group_by)
        data = plot.to_dict()['data']
        data = self.get_ordered_traces(data, group_by)
        for trace in data:
            # Removing nan that can not be converted to JSON
            trace['x'] = self.np2list(trace['x'])
            trace['y'] = self.np2list(trace['y'])
            trace['customdata'] = trace['customdata'].tolist()
        return data

    @classmethod
    def _get_reordered_traces(cls, traces: List[Dict], ordered_names: List[str]) -> List[Dict]:
        """
        Returns a new list of traces in the right order (defined in Categorical if available).
        Categorical can sometimes have unused categorizes (size equal or larger than traces)
        :param traces:
        :param ordered_names:
        :return:
        """
        if ordered_names and len(ordered_names) >= len(traces):
            lookup = {trace['name']: trace for trace in traces}
            ordered_traces = [lookup[name] for name in ordered_names if name in lookup]
            if len(ordered_traces) == len(traces):
                return ordered_traces
        return traces

    def _get_ordered_categories(self, group_by: str) -> List[str] | None:
        """
        Returns the categorical order for that obs field.
        :param group_by:
        :return:
        """
        if group_by in self.data.obs and self.data.obs[group_by].dtype.name == "category":
            return self.data.obs[group_by].cat.categories.astype(str).tolist()

    def get_ordered_traces(self, traces: List[Dict], group_by: str) -> List[Dict]:
        """
        Reorders traces based on user specified order in Categorical.
        :param traces:
        :param group_by:
        :return:
        """
        try:
            ordered_names = self._get_ordered_categories(group_by)
            return self._get_reordered_traces(traces, ordered_names)
        except Exception as e:
            logger.warning(f'Could not reorder traces: {e}')
            return traces

    def _get_default_info(self, field: str, na='n/a') -> str:
        """
        Sometimes a dataset might be missing some columns for info values. For example, it might not have a 'condition'
        column. In that, case, for all cells we will return the default value for 'condition', stored in
        'infoDefault' field of meta.

        :param field: Name of the info field, such as 'condition', 'sample', or 'cellType'
        :param na:
        :return:
        """
        return self.meta.get('infoDefault', {}).get(field) or na

    def get_info(self, field: str) -> pd.Series:
        if field in self.meta['info']:
            data_field = self.meta['info'][field]
            if data_field in self.data.obs:
                return self.data.obs[data_field].astype('str')
            else:
                raise UserException(f'No {data_field} data in dataset.')
        else:
            default_info = self._get_default_info(field)
            return pd.Series([default_info] * self.data.obs.index.size)

    def get_cell_type_field(self) -> str:
        """
        Returns the obs field containing cell types.
        :return:
        """
        return self.meta.get('info', {}).get('cellType')

    def generate_cell_ids(self) -> ndarray:
        """
        Generate cell ids, corresponding to indexes of cells in the matrix.
        The cell ids are independent of the group/trace.
        Cell ids start at 0. shape[0] gives the number of rows, corresponding to cells.
        @todo Use official cell ids, such as in orig.ident
        :return:
        """
        return np.arange(self.data.shape[0])

    def get_tooltip_data(self) -> List[ndarray | Series]:
        return [
            self.generate_cell_ids(),  # Cell Id
            self.get_info('cellType'),  # Cell Type
            self.get_info('condition'),  # Condition
            self.get_info('sample')  # Sample ID
        ]

    def get_default_group_by(self) -> str:
        """
        Returns the default column used for assigning colors to cells in the dataset.
        Meta dict contains a list of columns that can be used to group cells into traces and assign them colors.
        The first item in the list is considered default grouping. For example, 'cellType' will color
        (and group intro traces) cells based on their cell type column.
        :return:
        """
        if not self.meta['groups']:
            raise UserException(
                'Groups declaration in the metadata for the dataset is empty. Please declare at least one group.')

        return self.meta['groups'][0]['key']

    def get_default_embedding_key(self) -> Union[str, None]:
        """
        Returns the default embeddings matrix name. ex: X_umap
        If the field embeddings is available in meta and is not empty, returns key for first item
        Otherwise, returns last embedding
        :return:
        """
        if self.meta.get('embeddings'):
            return self.meta['embeddings'][0]['key']
        return list(self.data.obsm.keys())[-1]

    def as_json(self, embedding_name: str = None, group_by: str = None) -> Dict:
        """
        Returns a JSON used by the frontend to display a scatter plot of the cells.

        :param group_by: Used to color cells.
        :param embedding_name: matrix containing 2D coordinates of cells.
        :return: Plot object
        """
        if not embedding_name:
            embedding_name = self.get_default_embedding_key()

        if not group_by:
            group_by = self.get_default_group_by()

        data = self.as_plot_data(embedding_name, group_by)

        genes = self.get_genes()

        return {
            'id': self.dataset_id,
            'data': data,
            'meta': self.meta,
            'embedding': embedding_name,
            'group': group_by,
            'genes': genes,
            'interactions': self.has_interactions()
        }

    def get_umap(self, embedding_name: str, group_by: str) -> Dict:
        """
        Returns information necessary to draw a UMAP plot.
        Supports partial embeddings (missing cells have nan/null x/y coordinates.)
        :param group_by: Column with categorical or numeric values for each cell.
        :param embedding_name: Matrix with 2D coordinates of the cells (x,y)
        :return:
        """
        logger.info(f'Building UMAP from embedding {embedding_name}, grouped by {group_by}')
        x, y = self.get_embedding(embedding_name=embedding_name)
        colors = self.get_group_by(group_by)
        return {
            'x': self.np2list(x),
            'y': self.np2list(y),
            'colors': [color if color != 'nan' else None for color in colors]
        }

    def get_embedding(self, embedding_name: str = None) -> Tuple[ndarray[float], ndarray[float]]:
        """
        Returns the PCA embedding coordinates as a tuple of X, Y
        :param embedding_name:
        :return:
        """
        if len(self.data.obsm) == 0:
            raise DatasetException(f'The dataset {self.dataset_id} has no Embeddings.')

        if not embedding_name:
            embedding_name = list(self.data.obsm.keys())[-1]  # Last embedding

        if embedding_name not in self.data.obsm:
            raise DatasetException(f'Embedding {embedding_name} not available.')

        embedding = self.data.obsm[embedding_name]
        x = embedding[:, 0]
        y = embedding[:, 1]

        return x, y

    def get_group_by(self, group_by: str = None) -> ndarray[str]:
        """
        Returns the column used to assign cells/points into traces (colors).
        :param group_by:
        :return:
        """
        if not group_by:
            raise DatasetException('Group By column must be provided')

        if group_by not in self.data.obs:
            raise DatasetException(f'Group by column {group_by} does not exist for dataset.')

        return self.data.obs[group_by].astype('str')

    @classmethod
    def np2list(cls, arr: ndarray[float]) -> List[float]:
        arr_list = arr.tolist()
        nan_mask = np.isnan(arr)
        return [None if is_nan else value for is_nan, value in zip(nan_mask, arr_list)]

    def get_cell_idx(self, cell_id: str) -> int:
        if cell_id not in self.data.obs_names:
            raise UserException(f'Cell Id not found in dataset: {cell_id}')
        return self.data.obs_names.get_loc(cell_id)

    def get_gene_idx(self, gene_name: str) -> int:
        if gene_name not in self.data.var_names:
            raise UserException(f'Gene not found in dataset: {gene_name}')
        return self.data.var_names.get_loc(gene_name)

    def get_expression(self, cell_id: str, gene_name: str) -> float:
        return self.data[cell_id, gene_name].X.toarray()[0, 0]

    def get_gene_expression(self, gene_name: str) -> Dict[str, float]:
        gene_expression_column = self.data[:, gene_name].X.toarray().squeeze()
        return dict(zip(self.data.obs_names, gene_expression_column))

    def get_gene_expressions(self, gene_name: str) -> List[float]:
        """
        Returns a list of expression values for a gene.
        :param gene_name:
        :return:
        """
        gene_idx = self.get_gene_idx(gene_name)
        return self.data.X[:, gene_idx].toarray().flatten().tolist()

    def get_cell_expression(self, cell_id: str) -> Dict[str, float]:
        cell_expression_row = self.data[cell_id, :].X.toarray().squeeze()
        return dict(zip(self.data.var_names, cell_expression_row))

    def get_subset_for_genes(self, genes: List[str]) -> AnnData:
        return self.data[:, self.data.var_names.isin(genes)]

    def get_max_expression(self) -> float:
        """
        Returns the max expression value for all genes across all cells.
        Uses precomputed value if available.
        :return:
        """
        if 'max_expression' in self.data.uns:
            return self.data.uns['max_expression']
        logger.warning('Precomputed max expression not available. Computing in real time, which might be expensive.')
        return self.data.X.to_memory().max()

    def get_gene_expression_matrix(self, gene: str) -> ArrayView:
        return self.data[:, gene].X

    def get_max_gene_expression(self, gene: str) -> float:
        return float(self.get_gene_expression_matrix(gene).max())

    def get_normalized_gene_expression(self, gene: str, max_expression: float) -> List[float]:
        """
        If max expression in 0, it is replaced by 1 to avoid division by 0.
        :param gene:
        :param max_expression:
        :return:
        """
        return (self.get_gene_expression_matrix(gene) / (max_expression or 1)).toarray().flatten().tolist()

    @classmethod
    def filter_hvg(cls, adata: AnnData, hvg_threshold, receipt: Dict = None) -> AnnData:
        """
        @warning Some very downregulated genes can be lost.
        Better not to use this filtering, but if have to, need to have at least 4K genes left.
        :param receipt:
        :param adata:
        :param hvg_threshold:
        :return:
        """
        if adata.shape[1] > hvg_threshold:
            logger.info(f'Number of genes more than limit of {hvg_threshold}. Filtering by HVGs...')
            hvg_info = sc.pp.highly_variable_genes(adata, inplace=False)
            subset = adata[:, hvg_info['highly_variable']]
            Receipt.add_hvg_step(receipt)
            logger.info(f'After filtering by HVG, number of genes went from {adata.shape[1]} to {subset.shape[1]}')
            return subset

    @classmethod
    def subsample(cls, adata: AnnData, max_cells=1000, receipt: Dict = None):
        """
        The older, now deprecated subsampling method. Uses scanpy subsample.
        :param adata:
        :param max_cells:
        :param receipt:
        :return:
        """
        if adata.shape[0] > max_cells:
            logger.info(f'Number of cells more than limit of {max_cells}. Subsampling...')
            current_x_shape = adata.shape
            sc.pp.subsample(adata, n_obs=max_cells)
            Receipt.add_scanpy_subsample_step(receipt, n_obs=max_cells)
            logger.info(f'After subsampling, number of cells went from {current_x_shape[0]} to {adata.shape[0]}')

    def subsample_uniform(self, group_name: str, cell_ids: List[str], threshold: int = 1000,
                          receipt: Dict = None) -> List[str]:
        """
        Subsamples cell ids into uniform cluster representation.
        Cell Types are used for clusters.
        Final number of samples per cluster doesn't depend on original cluster sizes.
        An alternative would be a proportional split.

        :param group_name: Name of the group, used for receipt and logging.
        :param cell_ids: List of cell ids.
        :param threshold: Max number of cells allows in this group.
        :param receipt: [optional] Tracks data processing steps.
        :return: A list of cell ids (NOT indexes)
        """
        if len(cell_ids) <= threshold:
            return cell_ids

        # We group all cell indexes by cell type.
        cluster_field = self.get_cell_type_field()
        if not cluster_field:
            logger.warning(f'Dataset has no cell type field declared. Subsampling whole group.')
            return random.sample(cell_ids, threshold)

        clusters = self.get_clustered_cell_ids(cell_ids, cluster_field)

        # We want all cell ids to be split evenly by cell type. Each cell type will be limited by same limit.
        max_cells_per_cluster = int(math.floor(threshold / len(clusters)))
        original_clusters = clusters.copy()

        # Cutting every cluster to the limit if necessary. Using reproducible random.
        for cluster_name, cluster_cell_ids in clusters.items():
            if len(cluster_cell_ids) > max_cells_per_cluster:
                clusters[cluster_name] = random.sample(cluster_cell_ids, max_cells_per_cluster)

        # Into a flat list
        sampled_cell_indexes = [cell for cells in clusters.values() for cell in cells]

        # Adding to receipt
        Receipt.add_custom_subsample_step(receipt, group_name=group_name, group_max_cells=threshold,
                                          before_cells=cell_ids, after_cells=sampled_cell_indexes,
                                          before_clusters=original_clusters, after_clusters=clusters,
                                          max_cells_per_cluster=max_cells_per_cluster)

        return sampled_cell_indexes

    def get_clustered_cell_ids(self, cell_ids: List[str], cluster_col: str) -> Dict[str, List[str]]:
        """
        Returns cell indexes grouped by cluster
        :return:
        """
        # noinspection PyTypeChecker
        return self.data[cell_ids].obs.groupby(cluster_col, observed=True).apply(
            lambda df: df.index.tolist()).to_dict()

    def get_cell_ids(self, cell_indexes: List[int]) -> List[str]:
        return self.data.obs.index[cell_indexes].tolist()
