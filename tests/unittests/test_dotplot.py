from unittest import TestCase

from src.cellar3.analysis.differential.plots.dot import DotPlot
from src.cellar3.analysis.differential.subset import ExpressionMatrix
from tests.unittests.mock_dataset import MockDataset


class DotPlotTests(TestCase):

    def test_as_plot_data(self):
        target = [0, 2]
        background = [3]
        genes = ['Gene2']
        dataset = MockDataset('example1')
        adata = ExpressionMatrix.extract(dataset=dataset, target=target, background=background, meta={},
                                         min_cells_group=1)
        plot = DotPlot()
        plot_data = plot.as_plot_data(adata, genes=genes, meta={})
        self.assertListEqual(plot_data, [
            {'genes': ['Gene2', 'Gene2'],
             'groups': ['target', 'background'],
             'means': [1.0, 0.05],
             'percents': [1.0, 1.0]}
        ])
