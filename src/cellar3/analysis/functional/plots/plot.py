import abc
from typing import Dict, List

from pandas import DataFrame


class AbstractEnrichmentPlot:
    """
    A plot displaying data for specific enrichment categories.
    """

    @abc.abstractmethod
    def as_plot_data(self, df: DataFrame, meta: Dict[str, any]) -> List[Dict]:
        pass
