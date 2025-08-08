import logging
import warnings
from typing import List, Optional, Literal, Dict

import numpy as np
import scanpy as sc
from anndata import AnnData
from pandas import DataFrame

from src.cellar3.receipt import Receipt
from src.cellar3.tools import timeit

logger = logging.getLogger('cellar.analysis.expression.diff')

Method = Optional[Literal['logreg', 't-test', 'wilcoxon', 't-test_overestim_var']]
Category = Literal['non-significant', 'downregulated', 'upregulated']


class DifferentialExpression:
    """
    https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html

    Possible alternative for DE is the library diffxpy
    """

    CATEGORY_INSIGNIFICANT = 'non-significant'
    CATEGORY_DOWN_REGULATED = 'downregulated'
    CATEGORY_UP_REGULATED = 'upregulated'

    def __init__(self, method: Method = 'wilcoxon', pval_thr=0.05,
                 fold_change_threshold=0.25):
        """

        :param method: Using 'wilcoxon' by default to be similar to Seurat FindMarkers
        :param pval_thr:
        :param fold_change_threshold: Using 0.25 to have same default as FindMarkers
        """

        self.fold_change_thr = fold_change_threshold
        self.pval_thr = pval_thr
        self.method = method

    def assign_category(self, df: DataFrame):
        df['category'] = self.CATEGORY_INSIGNIFICANT
        df.loc[(df['pvals_adj'] < self.pval_thr) & (
                df['log2_fold_change'] > self.fold_change_thr), 'category'] = self.CATEGORY_UP_REGULATED
        df.loc[(df['pvals_adj'] < self.pval_thr) & (
                df['log2_fold_change'] < -self.fold_change_thr), 'category'] = self.CATEGORY_DOWN_REGULATED

    @classmethod
    def clean_data(cls, df: DataFrame):
        df['log2_fold_change'].replace([np.inf, -np.inf], None, inplace=True)
        df['log2_fold_change'].replace([np.nan], 0, inplace=True)  # no change or missing data
        df['minus_log10_pvals'].replace([np.inf, -np.inf, np.nan], None, inplace=True)

    @classmethod
    def _to_data_frame(cls, dge) -> DataFrame:
        df = DataFrame({
            'gene_names': dge['names']['target'],
            'log2_fold_change': dge['logfoldchanges']['target'],
            'pvals': dge['pvals']['target'],
            'pvals_adj': dge['pvals_adj']['target'],
        })

        df['minus_log10_pvals'] = -np.log10(df['pvals_adj'])

        return df

    @staticmethod
    def get_genes(df: DataFrame, category: str) -> List[str]:
        return df[df['category'] == category]['gene_names'].tolist()

    @staticmethod
    def get_up_genes(df: DataFrame) -> List[str]:
        return DifferentialExpression.get_genes(df, DifferentialExpression.CATEGORY_UP_REGULATED)

    @staticmethod
    def get_down_genes(df: DataFrame) -> List[str]:
        return DifferentialExpression.get_genes(df, DifferentialExpression.CATEGORY_DOWN_REGULATED)

    @timeit
    def process(self, adata: AnnData, target: List[int], background: List[int], meta: Dict[str, any]) -> DataFrame:
        """
        @todo Pass method as parameter
        :param adata: AnnData object that as already assigned groups and filtered
        :param target: group of interest, where we want to see up-regulated / down-regulated genes.
        :param background: Group to which we compare genes from target group.
        :param meta: Information on how the dataset was processed to get differential.
        :return:
        """
        logger.info(f'Received target and background. Target: {len(target)}. Background: {len(background)}')
        Receipt.add_dge_step(meta, method=self.method)
        Receipt.assign_category(meta, params={'p-val threshold': self.pval_thr,
                                              'fold change threshold': self.fold_change_thr})

        if adata:
            # scanpy does not allow to return the results without modifying the original adata object.
            # because the adata we get from extract is a view, we need to explicitely convert it to editable object.
            warnings.filterwarnings('ignore', category=FutureWarning, module='numpy')
            adata = adata.copy()
            sc.tl.rank_genes_groups(adata, groupby='group', groups=['target'], reference='background',
                                    method=self.method)
            dge = adata.uns['rank_genes_groups']  # differential is a Dict containing pandas objects
            df = self._to_data_frame(dge)
            self.assign_category(df)  # Assigning Category
            self.clean_data(df)  # Removing invalid numbers
            return df
