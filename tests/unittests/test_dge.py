from unittest import TestCase

from src.cellar3.analysis.differential.dge import DifferentialExpression
from src.cellar3.analysis.differential.subset import ExpressionMatrix
from tests.unittests.mock_dataset import MockDataset


class DifferentialGeneExpressionTests(TestCase):

    def test_process(self):
        target = [0, 2]
        background = [3, 4, 5]
        dataset = MockDataset('example2')
        adata = ExpressionMatrix.extract(dataset=dataset, target=target, background=background, meta={},
                                         min_cells_group=1, min_expression=0)
        dge = DifferentialExpression().process(adata, target=target, background=background, meta={})
        self.assertListEqual(dge['gene_names'].tolist(), ['Gene3', 'Gene2'])
        self.assertListEqual(dge['log2_fold_change'].tolist(), [-1.6872150897979736, 0.0])
        self.assertListEqual(dge['minus_log10_pvals'].tolist(), [-0.0, -0.0])
        self.assertListEqual(dge['category'].tolist(), ['non-significant', 'non-significant'])
