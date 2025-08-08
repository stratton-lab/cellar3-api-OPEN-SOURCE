import logging
from typing import Dict, Any

from pandas import DataFrame

from src.cellar3.analysis.functional.plots.plot import AbstractEnrichmentPlot
from src.cellar3.receipt import Receipt

logger = logging.getLogger('cellar.enrichment.plot.bar')


class EnrichmentBarPlot(AbstractEnrichmentPlot):

    def __init__(self, top_terms=10):
        self.top_terms = top_terms

    def as_plot_data(self, df: DataFrame, meta: Dict[str, any]) -> Any:
        Receipt.bar_plot(meta, params={'top_terms': self.top_terms, 'sort': 'minus_log10_pvals DESC'})
        return df[['Term', 'Combined Score', 'P-value', 'minus_log10_pvals']] \
            .sort_values(by='minus_log10_pvals', ascending=False) \
            .head(self.top_terms) \
            .to_dict(orient='records')
