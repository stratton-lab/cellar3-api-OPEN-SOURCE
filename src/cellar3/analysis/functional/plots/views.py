import logging

import numpy as np
from rest_framework.response import Response
from rest_framework.views import APIView

from src.cellar3.analysis.differential.cached import get_cached_dge
from src.cellar3.analysis.functional.enrich import EnrichmentAnalysis, Enrichir
from src.cellar3.analysis.functional.plots.bar import EnrichmentBarPlot
from src.cellar3.analysis.functional.plots.network import EnrichmentNetworkPlot
from src.cellar3.analysis.functional.plots.plot import AbstractEnrichmentPlot
from src.cellar3.exceptions import UserException
from src.cellar3.receipt import Receipt

logger = logging.getLogger('cellar.analysis.functional.plots')


class AbstractEnrichmentView(APIView):
    analyzer: EnrichmentAnalysis
    plot: AbstractEnrichmentPlot

    def get(self, request, dataset_id):
        return Response({'status': 'This plot is only available through POST'})

    def post(self, request, dataset_id):
        try:
            meta = Receipt.create()

            # Input
            json_data = request.data
            target = json_data.get('target')
            background = json_data.get('background')
            gene_set = json_data.get('gene_set')
            diff_p_val_thr = 0.05  # P-value Threshold
            fc_thr = 0.25  # Fold Change Threshold

            dataset, meta, dge = get_cached_dge(dataset_id, target, background, diff_p_val_thr, fc_thr, meta)
            enrichment = self.analyzer.process(gene_set=gene_set, dge=dge, organism=dataset.meta['species'],
                                               meta=meta)
            enrichment['minus_log10_pvals'] = -np.log10(enrichment['Adjusted P-value'])
            traces = self.plot.as_plot_data(df=enrichment, meta=meta)
            term2genes = Enrichir.get_term_genes(dge, enrichment)

            return Response({'data': traces, 'meta': meta, 'term2genes': term2genes})
        except Exception as e:
            logger.exception(e)
            raise UserException(f'Could not process Enrichment: {e}')


class BarPlotView(AbstractEnrichmentView):
    analyzer = Enrichir()
    plot = EnrichmentBarPlot()


class NetworkPlotView(AbstractEnrichmentView):
    analyzer = Enrichir()  # Prerank() Slower and current has bug
    plot = EnrichmentNetworkPlot()
