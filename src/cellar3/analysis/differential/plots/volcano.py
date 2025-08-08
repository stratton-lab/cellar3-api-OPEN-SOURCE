import logging
from typing import Dict, List

from pandas import DataFrame

from src.cellar3.analysis.differential.dge import DifferentialExpression
from src.cellar3.analysis.differential.plots.plot import DifferentialGeneExpressionPlot

logger = logging.getLogger('cellar.plot.volcano')


class VolcanoPlot(DifferentialGeneExpressionPlot):

    @classmethod
    def _get_trace_data(cls, df: DataFrame, category: str) -> Dict:
        return {
            'category': category,
            'genes': df['gene_names'].tolist(),
            'logFC': df['log2_fold_change'].tolist(),
            'minusLogP': df['minus_log10_pvals'].tolist(),
        }

    def as_plot_data(self, dge: DataFrame, meta: Dict[str, any]) -> List[Dict]:
        """
        Generates data used by front-end to draw a Volcano plot.

        :param meta:
        :param dge: Differential Gene Expression data calculated from the main dataset.
        :return: data used by front-end to draw a Volcano plot.
        """
        return [self._get_trace_data(dge[dge['category'] == category], category) for category in
                [
                    DifferentialExpression.CATEGORY_INSIGNIFICANT,
                    DifferentialExpression.CATEGORY_DOWN_REGULATED,
                    DifferentialExpression.CATEGORY_UP_REGULATED
                ]]
