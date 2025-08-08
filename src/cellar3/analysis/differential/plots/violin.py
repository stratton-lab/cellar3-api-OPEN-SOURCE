import logging
from typing import List, Dict

import numpy as np
from anndata import AnnData

from src.cellar3.analysis.differential.plots.plot import AbstractGenePlot
from src.cellar3.analysis.differential.subset import GROUP_COL, GROUP_BACKGROUND, GROUP_TARGET
from src.cellar3.exceptions import UserException

logger = logging.getLogger('cellar.analysis.expression.violin')


class ViolinPlot(AbstractGenePlot):

    def as_plot_data(self, adata: AnnData, genes: List[str], meta: Dict[str, any]) -> List[Dict]:
        if len(genes) == 0:
            raise UserException('No genes selected.')

        if adata:
            traces = []
            group_col = GROUP_COL
            group_data = adata.obs[group_col]

            for group in [GROUP_BACKGROUND, GROUP_TARGET]:
                if group not in group_data.unique():
                    logger.info(f"Group {group} not found in adata.obs. Skipping.")
                    continue

                side = 'positive' if group == 'target' else 'negative'
                group_indices = group_data[group_data == group].index
                integer_group_indices = [adata.obs_names.get_loc(id) for id in group_indices]
                logger.info(f'Group {group} has {len(integer_group_indices)} sample ids.')
                trace = {'name': group, 'side': side, 'x': [], 'y': []}
                for gene in genes:
                    matches = np.where(adata.var_names == gene)[0]
                    if not len(matches):
                        raise UserException(f'Gene {gene} not found in subset.')
                    gene_index = matches[0]
                    gene_data = adata.X[integer_group_indices, gene_index].toarray().flatten()
                    logger.info(
                        f'Group {group} Gene {gene} Expression: mean({np.mean(gene_data)}, max({np.max(gene_data)})')
                    trace['x'].extend([gene] * len(gene_data))
                    trace['y'].extend(gene_data.tolist())
                traces.append(trace)

            return traces
