from unittest import TestCase

from src.cellar3.analysis.differential.subset import ExpressionMatrix, GROUP_UNASSIGNED
from tests.unittests.mock_dataset import MockDataset
from tests.unittests.tools import get_group


class ExpressionMatrixTests(TestCase):
    def test_check_parameters(self):
        target = [0, 2]
        background = [3]
        ExpressionMatrix._check_parameters(target=target, background=background)

    def test_assign_groups(self):
        adata = MockDataset('example1').data
        self.assertEqual(adata.shape, (6, 3))
        target = ['Cell1_T', 'Cell3_T']
        background = ['Cell4_B']
        ExpressionMatrix._assign_groups(adata, target=target, background=background)
        self.assertEqual(get_group(adata, 'Cell1_T'), 'target')
        self.assertEqual(get_group(adata, 'Cell2'), 'unassigned')
        self.assertEqual(get_group(adata, 'Cell3_T'), 'target')
        self.assertEqual(get_group(adata, 'Cell4_B'), 'background')
        self.assertEqual(get_group(adata, 'Cell5'), 'unassigned')
        self.assertEqual(get_group(adata, 'Cell6'), 'unassigned')

    def test_check_subset(self):
        target = ['Cell1_T', 'Cell3_T']
        background = ['Cell4_B']
        adata = MockDataset('example1').data
        ExpressionMatrix._assign_groups(adata, target=target, background=background)
        subset = adata[adata.obs['group'] != GROUP_UNASSIGNED].to_memory()
        ExpressionMatrix._check_subset(adata, subset, target, background)
        self.assertEqual(subset.shape, (3, 3))
        self.assertListEqual(subset.obs.index.tolist(), ['Cell1_T', 'Cell3_T', 'Cell4_B'])

    def test_filter_groups(self):
        target = ['Cell1_T', 'Cell3_T']
        background = ['Cell4_B']
        adata = MockDataset('example1').data
        ExpressionMatrix._assign_groups(adata, target=target, background=background)
        filtered_by_cells = adata[adata.obs['group'] != GROUP_UNASSIGNED].to_memory()

        filtered_by_genes = ExpressionMatrix._filter_groups(filtered_by_cells, min_cells_group=1,
                                                            min_expression=0)
        self.assertEqual(filtered_by_genes.shape, (3, 1))
        self.assertListEqual(filtered_by_genes.obs.index.tolist(), ['Cell1_T', 'Cell3_T', 'Cell4_B'])
        self.assertListEqual(filtered_by_genes.var.index.tolist(), ['Gene2'])

        filtered_by_genes = ExpressionMatrix._filter_groups(filtered_by_cells, min_cells_group=2,
                                                            min_expression=0)
        self.assertEqual(filtered_by_genes.shape, (3, 0))
        self.assertListEqual(filtered_by_genes.obs.index.tolist(), ['Cell1_T', 'Cell3_T', 'Cell4_B'])
        self.assertListEqual(filtered_by_genes.var.index.tolist(), [])

    def test_preprocess(self):
        target = [0, 2]
        background = [3]
        dataset = MockDataset('example1')
        preprocessed_adata = ExpressionMatrix.extract(
            dataset, target=target, background=background, meta={},
            min_cells_group=1, min_expression=0)
        self.assertEqual(preprocessed_adata.shape, (3, 1))
        self.assertListEqual(preprocessed_adata.obs.index.tolist(), ['Cell1_T', 'Cell3_T', 'Cell4_B'])
        self.assertListEqual(preprocessed_adata.var.index.tolist(), ['Gene2'])
