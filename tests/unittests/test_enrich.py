from unittest import TestCase

from src.cellar3.analysis.differential.dge import DifferentialExpression
from src.cellar3.analysis.differential.subset import ExpressionMatrix
from src.cellar3.analysis.functional.enrich import Enrichir
from tests.unittests.mock_dataset import MockDataset


class FunctionalEnrichmentTests(TestCase):

    def DEACTIVATED_test_process(self):
        target = list(range(0, 51))
        background = list(range(51, 100))
        dataset = MockDataset('example3')
        adata = ExpressionMatrix.extract(dataset=dataset, target=target, background=background, meta={},
                                         min_cells_group=1, min_expression=0)
        dge = DifferentialExpression(pval_thr=0.05, fold_change_threshold=0.025) \
            .process(adata, target=target, background=background, meta={})
        enrichment = Enrichir().process(gene_set='GO_Biological_Process_2023', dge=dge, organism='Human', meta={})
        self.assertEqual(enrichment['Term'][0], 'Positive Regulation Of Transmembrane Transport (GO:0034764)')
