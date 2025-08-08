from unittest import TestCase

from src.cellar3.analysis.differential.dge import DifferentialExpression
from src.cellar3.analysis.differential.plots.volcano import VolcanoPlot
from src.cellar3.analysis.differential.subset import ExpressionMatrix
from tests.unittests.mock_dataset import MockDataset


class VolcanoPlotTests(TestCase):

    def test_as_plot_data(self):
        target = [0, 2]
        background = [3, 4, 5]
        dataset = MockDataset('example2')
        adata = ExpressionMatrix.extract(dataset=dataset, target=target, background=background, meta={},
                                         min_cells_group=1, min_expression=0)
        dge = DifferentialExpression().process(adata, target=target, background=background, meta={})
        traces = VolcanoPlot().as_plot_data(dge=dge, meta={})
        self.assertListEqual(traces, [
            {'category': 'non-significant',
             'genes': ['Gene3', 'Gene2'],
             'logFC': [-1.6872150897979736, 0.0],
             'minusLogP': [-0.0, -0.0]},
            {'category': 'downregulated', 'genes': [], 'logFC': [], 'minusLogP': []},
            {'category': 'upregulated', 'genes': [], 'logFC': [], 'minusLogP': []}
        ])
