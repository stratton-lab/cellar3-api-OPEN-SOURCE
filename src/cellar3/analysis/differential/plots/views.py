import logging

from rest_framework.response import Response
from rest_framework.views import APIView

from src.cellar3.analysis.differential.dge import DifferentialExpression
from src.cellar3.analysis.differential.plots.dot import DotPlot
from src.cellar3.analysis.differential.plots.plot import AbstractGenePlot
from src.cellar3.analysis.differential.plots.umap import UMAPPlot
from src.cellar3.analysis.differential.plots.violin import ViolinPlot
from src.cellar3.analysis.differential.plots.volcano import VolcanoPlot
from src.cellar3.analysis.differential.subset import ExpressionMatrix
from src.cellar3.datasets.dataset import Dataset
from src.cellar3.exceptions import UserException
from src.cellar3.receipt import Receipt

logger = logging.getLogger('cellar.analysis.plots')


class SimplePlotView(APIView):
    plot_generator: AbstractGenePlot

    def get(self, request, dataset_id):
        return Response({'status': 'This plot can only be generated through POST'})

    def post(self, request, dataset_id):
        try:
            meta = Receipt.create()
            json_data = request.data
            genes = json_data.get('genes')
            target = json_data.get('target')
            background = json_data.get('background')
            dataset = Dataset(dataset_id)
            # We set min_cells_group=0 because we do not want to filter by genes for plots.
            adata = ExpressionMatrix.extract(dataset=dataset, target=target, background=background, meta=meta,
                                             min_cells_group=0, genes=genes)
            traces = self.plot_generator.as_plot_data(adata=adata, genes=genes, meta=meta)
            return Response({'traces': traces, 'receipt': meta})
        except Exception as e:
            logger.exception(e)
            raise UserException(f'Could not generate plot: {e}')


class ViolinPlotView(SimplePlotView):
    plot_generator = ViolinPlot()


class DotPlotView(SimplePlotView):
    plot_generator = DotPlot()


class UmapPlotView(APIView):
    def get(self, request, dataset_id):
        return Response({'status': 'This plot can only be generated through POST'})

    def post(self, request, dataset_id):
        try:
            meta = Receipt.create()
            json_data = request.data
            embedding_name = json_data.get('embedding', None)
            genes = json_data.get('genes')
            normalize_to_gene = True
            dataset = Dataset(dataset_id)
            data = UMAPPlot().as_plot_data(dataset=dataset, embedding_name=embedding_name, genes=genes,
                                           normalize_to_gene=normalize_to_gene, meta=meta)
            return Response(data)
        except Exception as e:
            logger.exception(e)
            raise UserException(f'Could not generate plot: {e}')


class VolcanoPlotView(APIView):

    def get(self, request, dataset_id):
        return Response({'status': 'This plot is only available through POST'})

    def post(self, request, dataset_id):
        try:
            receipt = Receipt.create()
            json_data = request.data
            target = json_data.get('target')
            background = json_data.get('background')
            dataset = Dataset(dataset_id)
            Receipt.add_dataset_load_step(receipt, dataset.meta)
            adata = ExpressionMatrix.extract(dataset, target=target, background=background, meta=receipt)
            dge = DifferentialExpression().process(adata, target=target, background=background, meta=receipt)
            traces = VolcanoPlot().as_plot_data(dge=dge, meta=receipt)
            # significant_genes = SignificantGenes().extract(dge)
            return Response({'traces': traces, 'meta': receipt})
        except Exception as e:
            logger.exception(e)
            raise UserException(f'Could not process Differential Gene Expression: {e}')
