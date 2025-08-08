import logging
from typing import List, Dict

from anndata import AnnData

from src.cellar3.exceptions import UserException
from src.cellar3.receipt import Receipt
from src.cellar3.tools import timeit

logger = logging.getLogger('cellar.analysis.subset')

GROUP_COL = 'group'
GROUP_UNASSIGNED = 'unassigned'
GROUP_TARGET = 'target'
GROUP_BACKGROUND = 'background'


class ExpressionMatrix:
    HVG_FILTER_TRIGGER = 10000  # We use hgv filtering if number of genes is over the limit.
    SUBSAMPLE_TRIGGER = 1000  # Reduces the number of cells.

    @classmethod
    @timeit
    def extract(cls, dataset, target: List[int], background: List[int], meta: Dict[str, any],
                min_cells_group: int = 3, min_expression: float = 0, genes: List[str] = None) -> AnnData:
        """
        Assigning cells to groups, removing unassigned
        @note: Removed filtering min cells per gene as only works if dataset has raw counts
        @todo If len of target same as len of all cells in dataset and no background, skip filter by cells
        :param genes: Optional list of genes to prefilter the matrix.
        :param min_cells_group: Minimum number of cells in which a gene is expressed, necessary to keep that gene.
        :param min_expression: Minimum expression level of a gene in cell to consider it to be expressed.
        :param meta:
        :param dataset:
        :param target:
        :param background:
        :return:
        """
        adata = dataset.data

        if not adata.shape[0]:
            raise UserException("The initial dataset is empty.")

        cls._check_parameters(target, background)

        # Converting cell indexes to ids
        target = dataset.get_cell_ids(target)
        background = dataset.get_cell_ids(background)

        # Subsampling
        target = dataset.subsample_uniform('target', target, receipt=meta)
        background = dataset.subsample_uniform('background', background, receipt=meta)

        cls._assign_groups(adata, target, background)

        # Keep only cells from target & background
        subset = cls._load_subset_to_memory(adata, genes)
        # Makes sure we have cells left after filtering by target & background
        cls._check_subset(adata, subset, target, background)

        # Removes genes expressed in too few cells in either background or target
        # todo Only if min_cells_group > 0
        subset = cls._filter_groups(subset, min_cells_group, min_expression, receipt=meta)

        return subset

    @classmethod
    @timeit
    def _load_subset_to_memory(cls, adata: AnnData, genes: List[str]) -> AnnData:
        """
        Loads necessary expression sub-matrix into memory. This is the SLOWEST part.
        It is necessary because scanpy DGE needs all data into memory.
        :param adata: On disk expression data
        :return: Subset of adata fully loaded into RAM.
        """
        gene_filter = adata.var_names.isin(genes) if genes else slice(None)
        cell_filter = adata.obs['group'] != GROUP_UNASSIGNED
        subset = adata[cell_filter, gene_filter]
        return subset.to_memory()

    @classmethod
    @timeit
    def _filter_groups(cls, subset: AnnData, min_cells_group: int, min_expression: float,
                       receipt: Dict = None) -> AnnData:
        """
        Only keeps genes that are expressed in more that n cells in both target and background.
        :param subset:
        :param min_cells_group: Minimum number of cells in which a gene is expressed, necessary to keep that gene.
        :param min_expression: Minimum expression level of a gene in cell to consider it to be expressed.
        :return:
        """
        target = cls._genes_expressed_in_group(subset, GROUP_TARGET, min_cells_group, min_expression)
        background = cls._genes_expressed_in_group(subset, GROUP_BACKGROUND, min_cells_group, min_expression)
        filtered_subset = subset[:, (target & background)]

        logger.info(f'Only keeping genes expressed in at least {min_cells_group} cells in both target and background.')
        logger.info(f'After filtering, number of genes went from {subset.shape[1]} to {filtered_subset.shape[1]}')
        Receipt.add_gene_filter_step(receipt, min_cells_group)
        return filtered_subset

    @classmethod
    def _genes_expressed_in_group(cls, subset, group: str, min_cells_group: int, expression_threshold: float):
        return (subset[subset.obs[GROUP_COL] == group].X > expression_threshold).sum(axis=0) >= min_cells_group

    @classmethod
    def _check_parameters(cls, target: List[int], background: List[int]):
        if not target:
            raise UserException('No cells for Target')

        # if not background:
        #    raise UserException('No cells for Background')

        if set(target) & set(background):
            raise UserException("You can not have the same cells in both Target and Background.")

    @classmethod
    @timeit
    def _assign_groups(cls, adata, target: List[str], background: List[str]):
        """
        Assigns cells in adata to one of categories: unassigned, target, background.
        NOTE: We need to pass cell ids, NOT cell indexes
        :param adata:
        :param target: List of cell ids
        :param background: List of cell ids
        :return:
        """
        adata.obs['group'] = GROUP_UNASSIGNED
        adata.obs.loc[target, 'group'] = GROUP_TARGET
        adata.obs.loc[background, 'group'] = GROUP_BACKGROUND
        adata.obs['group'] = adata.obs['group'].astype('category')
        adata.strings_to_categoricals()

    @classmethod
    def _check_subset(cls, adata, subset, target, background):
        if subset.shape[0] == 0:
            raise UserException("Target and Background match no cells.")
        if subset.shape[0] != len(target) + len(background):
            raise UserException("Target and Background did not match all the cells.")

        logger.info(f'After loading filtered subset into memory, number of cells went from'
                    f' {adata.shape[0]} to {subset.shape[0]} and number of genes went '
                    f'from {adata.shape[1]} to {subset.shape[1]}')
