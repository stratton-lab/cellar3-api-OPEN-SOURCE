import io
import logging
from typing import Dict, List

import numpy as np
import scanpy as sc
from anndata import AnnData
from matplotlib import pyplot as plt

from src.cellar3.analysis.differential.plots.plot import AbstractGenePlot
from src.cellar3.analysis.differential.subset import GROUP_COL
from src.cellar3.exceptions import UserException

logger = logging.getLogger('cellar.analysis.expression.dot')


class DotPlot(AbstractGenePlot):

    def __init__(self, expression_threshold=0):
        self.expression_threshold = expression_threshold

    def as_plot_data(self, adata: AnnData, genes: List[str], meta: Dict[str, any]) -> List[Dict]:
        """
        Each dot represents cells for a gene in a group.
        - Size: percentage of cells expressing the feature in each cluster.
        - Color: average expression level

        # WARNING: Will not work if expression normalized and have negative values
        :param genes:
        :param adata:
        :param meta:
        :return:
        """
        if len(genes) == 0:
            raise UserException('No genes selected.')

        if adata:
            logger.info(f'Initial matrix has {adata.X.shape[0]} cells and {adata.X.shape[1]} genes.')
            gene_mask = adata.var_names.isin(genes)
            groups = adata.obs[GROUP_COL].unique()
            trace = {'groups': [], 'genes': [], 'means': [], 'percents': []}

            for group in groups:
                group_mask = adata.obs['group'] == group
                view = adata[group_mask, gene_mask]

                pct_cells_expr = np.array((view.X > self.expression_threshold).mean(axis=0)).flatten()
                avg_expr = np.array(view.X.mean(axis=0)).flatten()  # Average expression level
                for gene in genes:
                    if gene in view.var_names:
                        idx = list(view.var_names).index(gene)
                        trace['genes'].append(gene)
                        trace['groups'].append(group)
                        trace['means'].append(avg_expr[idx])
                        trace['percents'].append(pct_cells_expr[idx])

            return [trace]

    def as_image(self, adata: AnnData, genes: List[str], meta: Dict[str, any]) -> bytes:
        """
        Returns a binary image that can be served with `return HttpResponse(data, content_type='image/png')`
        :param adata:
        :param genes:
        :param meta:
        :return:
        """
        if len(genes) == 0:
            raise UserException('No genes selected.')

        if adata:
            res = sc.pl.dotplot(adata, genes, groupby="group",
                                expression_cutoff=self.expression_threshold, return_fig=True)
            fig = res.fig
            buffer = io.BytesIO()
            fig.savefig(buffer, format='png', dpi=300)
            plt.close(fig)
            buffer.seek(0)
            return buffer.getvalue()
