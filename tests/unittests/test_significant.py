from unittest import TestCase

from src.cellar3.analysis.differential.dge import DifferentialExpression
from src.cellar3.analysis.differential.significant import SignificantGenes
from src.cellar3.analysis.differential.subset import ExpressionMatrix
from tests.unittests.mock_dataset import MockDataset


class SignificantGenesTests(TestCase):

    def test_extract(self):
        target = [0, 2]
        background = [3, 4, 5]
        dataset = MockDataset('example2')
        adata = ExpressionMatrix.extract(dataset=dataset, target=target, background=background, meta={},
                                         min_cells_group=1, min_expression=0)
        dge = DifferentialExpression().process(adata, target=target, background=background, meta={})
        significant_genes = SignificantGenes().extract(dge)
        self.assertListEqual(significant_genes, [])
