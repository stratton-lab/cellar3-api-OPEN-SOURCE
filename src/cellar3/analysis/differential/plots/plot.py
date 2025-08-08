import abc
from typing import Dict, List

from anndata import AnnData
from pandas import DataFrame


class AbstractGenePlot:
    """
    A plot displaying data for specific genes.
    """

    @abc.abstractmethod
    def as_plot_data(self, adata: AnnData, genes: List[str], meta: Dict[str, any]) -> List[Dict]:
        pass

    def as_image(self, adata: AnnData, genes: List[str], meta: Dict[str, any]) -> bytes:
        pass


class DifferentialGeneExpressionPlot:
    def as_plot_data(self, dge: DataFrame, meta: Dict[str, any]) -> Dict:
        pass
