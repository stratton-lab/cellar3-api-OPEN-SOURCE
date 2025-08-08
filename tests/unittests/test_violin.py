from unittest import TestCase

from src.cellar3.analysis.differential.plots.violin import ViolinPlot
from src.cellar3.analysis.differential.subset import ExpressionMatrix
from tests.unittests.mock_dataset import MockDataset


class ViolinPlotTests(TestCase):

    def test_as_plot_data(self):
        target = [0, 2]
        background = [3]
        genes = ['Gene2']
        dataset = MockDataset('example1')
        adata = ExpressionMatrix.extract(dataset=dataset, target=target, background=background, meta={},
                                         min_cells_group=1)
        plot = ViolinPlot()
        traces = plot.as_plot_data(adata, genes=genes, meta={})
        self.assertListEqual(traces, [
            {'name': 'background',
             'side': 'negative',
             'x': ['Gene2'],
             'y': [0.05]},
            {'name': 'target',
             'side': 'positive',
             'x': ['Gene2', 'Gene2'],
             'y': [1.0, 1.0]}
        ])
