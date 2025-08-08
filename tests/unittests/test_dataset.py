from unittest import TestCase

from tests.unittests.mock_dataset import MockDataset


class DatasetsTests(TestCase):

    def test_get_genes_count(self):
        dataset = MockDataset('example2')
        self.assertEqual(dataset.get_genes_count(), 3)

    def test_get_cells_count(self):
        dataset = MockDataset('example2')
        self.assertEqual(dataset.get_cells_count(), 6)

    def test_get_available_embeddings(self):
        dataset = MockDataset('example2')
        self.assertEqual(dataset.get_available_embeddings(), ['X_mock_embedding'])

    def test_qc(self):
        # TODO Implement
        pass

    def test_as_plot_data(self):
        self.maxDiff = None
        dataset = MockDataset('example2', meta={'info': {}})
        self.assertListEqual(dataset.as_plot_data(embedding_name='X_mock_embedding', group_by='DATsubtype'), [
            {'customdata': [[0, 'n/a', 'n/a', 'n/a'], [1, 'n/a', 'n/a', 'n/a'], [2, 'n/a', 'n/a', 'n/a']],
             'hovertemplate': 'color=Neuron<br>x=%{x}<br>y=%{y}<extra></extra>',
             'legendgroup': 'Neuron', 'marker': {'color': '#636efa', 'symbol': 'circle'},
             'mode': 'markers', 'name': 'Neuron', 'orientation': 'v', 'showlegend': True,
             'type': 'scatter', 'x': [0.1, 0.3, 0.5], 'xaxis': 'x', 'y': [0.2, 0.4, 0.6],
             'yaxis': 'y'},
            {'customdata': [[3, 'n/a', 'n/a', 'n/a'], [5, 'n/a', 'n/a', 'n/a']],
             'hovertemplate': 'color=nan<br>x=%{x}<br>y=%{y}<extra></extra>',
             'legendgroup': 'nan', 'marker': {'color': '#EF553B', 'symbol': 'circle'}, 'mode': 'markers',
             'name': 'nan', 'orientation': 'v', 'showlegend': True, 'type': 'scatter',
             'x': [0.7, 1.2], 'xaxis': 'x', 'y': [0.8, 1.3], 'yaxis': 'y'},
            {'customdata': [[4, 'n/a', 'n/a', 'n/a']],
             'hovertemplate': 'color=Macrophage<br>x=%{x}<br>y=%{y}<extra></extra>',
             'legendgroup': 'Macrophage', 'marker': {'color': '#00cc96', 'symbol': 'circle'},
             'mode': 'markers', 'name': 'Macrophage', 'orientation': 'v', 'showlegend': True,
             'type': 'scatter', 'x': [0.9], 'xaxis': 'x', 'y': [1.1], 'yaxis': 'y'}
        ])

    def test__get_default_info(self):
        dataset = MockDataset('example2')
        dataset.meta = {"infoDefault": {
            "condition": "Healthy Adult",
            "sample": "Colon"
        }}
        self.assertEqual(dataset._get_default_info('condition'), 'Healthy Adult')
        self.assertEqual(dataset._get_default_info('sample'), 'Colon')
        self.assertEqual(dataset._get_default_info('test'), 'n/a')

    def test__get_info(self):
        dataset = MockDataset('example2')
        dataset.meta = {
            "info": {"cellType": "DATsubtype"},
            "infoDefault": {"condition": "Healthy Adult", "sample": "Colon"}
        }
        self.assertListEqual(dataset.get_info('cellType').tolist(),
                             ['Neuron', 'Neuron', 'Neuron', 'nan', 'Macrophage', 'nan'])
        self.assertListEqual(dataset.get_info('condition').tolist(),
                             ['Healthy Adult', 'Healthy Adult', 'Healthy Adult', 'Healthy Adult', 'Healthy Adult',
                              'Healthy Adult'])
        self.assertListEqual(dataset.get_info('test').tolist(), ['n/a', 'n/a', 'n/a', 'n/a', 'n/a', 'n/a'])

    def test__get_cell_id(self):
        dataset = MockDataset('example2')
        self.assertListEqual(list(dataset.generate_cell_ids()), [0, 1, 2, 3, 4, 5])

    def test_get_tooltip_data(self):
        dataset = MockDataset('example2')
        dataset.meta = {
            "info": {"cellType": "DATsubtype"},
            "infoDefault": {"condition": "Healthy Adult", "sample": "Colon"}
        }
        tooltip_data = dataset.get_tooltip_data()
        self.assertEqual(len(tooltip_data), 4)
        self.assertListEqual(list(dataset.get_tooltip_data()[0]), [0, 1, 2, 3, 4, 5])

    def test_get_cell_idx(self):
        dataset = MockDataset('example2')
        self.assertEqual(dataset.get_cell_idx('Cell2'), 1)

    def test_get_gene_idx(self):
        dataset = MockDataset('example2')
        self.assertEqual(dataset.get_gene_idx('Gene1'), 0)

    def test_get_expression(self):
        dataset = MockDataset('example2')
        self.assertEqual(dataset.get_expression('Cell2', 'Gene1'), 1.0)

    def test_get_gene_expression(self):
        dataset = MockDataset('example2')
        self.assertDictEqual(dataset.get_gene_expression('Gene1'), {
            'Cell1_T': 0.0, 'Cell2': 1.0, 'Cell3_T': 0.5, 'Cell4_B': 0.0, 'Cell5_B': 0.0, 'Cell6_B': 0.0
        })

    def test_get_cell_expression(self):
        dataset = MockDataset('example2')
        self.assertDictEqual(dataset.get_cell_expression('Cell2'), {'Gene1': 1.0, 'Gene2': 1.0, 'Gene3': 2.0})

    def test_get_subset_for_genes(self):
        dataset = MockDataset('example2')
        subset = dataset.get_subset_for_genes(['Gene1', 'Gene3'])
        self.assertEqual(subset.shape, (6, 2))
        self.assertEqual(subset.X.todense().squeeze().tolist(), [
            [0.0, 1.0],  # Cell1_T => Gene1: 0,  Gene3: 1
            [1.0, 2.0],  # Cell2 => Gene1: 1, Gene3: 2
            [0.5, 10.0],  # Cell3_T
            [0.0, 0.0],  # Cell4_B
            [0.0, 0.0],  # Cell5_B
            [0.0, 20.0]  # Cell6_B => Gene1: 0, Gene3: 20
        ])
        self.assertEqual(subset.X.todense().flatten().tolist(), [[
            0.0, 1.0,
            1.0, 2.0,
            0.5, 10.0,
            0.0, 0.0,
            0.0, 0.0,
            0.0, 20.0]
        ])

    def test_get_max_expression(self):
        dataset = MockDataset('example2')
        self.assertEqual(dataset.get_max_expression(), 20)

    def test_get_reordered_traces(self):
        dataset = MockDataset('example2')
        original_traces = [{'name': 'A'}, {'name': 'B'}, {'name': 'C'}]
        self.assertListEqual(dataset._get_reordered_traces(traces=original_traces, ordered_names=[]), original_traces)
        self.assertListEqual(dataset._get_reordered_traces(traces=original_traces, ordered_names=['B', 'C', 'A']),
                             [{'name': 'B'}, {'name': 'C'}, {'name': 'A'}])
        # If not enough info to reorder all, we don't reorder (fewer categories than actual values)
        self.assertListEqual(dataset._get_reordered_traces(traces=original_traces, ordered_names=['B', 'C']),
                             original_traces)
        # Some categories are not used
        self.assertListEqual(dataset._get_reordered_traces(traces=original_traces, ordered_names=['B', 'C', 'A', 'D']),
                             [{'name': 'B'}, {'name': 'C'}, {'name': 'A'}])
