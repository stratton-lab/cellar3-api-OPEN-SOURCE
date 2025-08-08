from unittest import TestCase

from src.cellar3.analysis.differential.plots.umap import UMAPPlot
from tests.unittests.mock_dataset import MockDataset


class UMAPPlotTests(TestCase):

    def test_as_plot_data(self):
        dataset = MockDataset('example2')
        plot = UMAPPlot()
        genes = ['Gene2']
        plot_data = plot.as_plot_data(dataset, None, genes=genes, meta={})
        self.assertListEqual(plot_data['genes'][0]['expression'], [0.05, 0.05, 0.05, 0.05, 0.05, 0.05])
        self.assertEqual(plot_data['genes'][0]['gene'], 'Gene2')
        self.assertListEqual(list(plot_data['x']), [0.1, 0.3, 0.5, 0.7, 0.9, 1.2])
        self.assertListEqual(list(plot_data['y']), [0.2, 0.4, 0.6, 0.8, 1.1, 1.3])
