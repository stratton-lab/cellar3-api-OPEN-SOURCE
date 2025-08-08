import logging
from typing import Dict

from pandas import DataFrame

from src.cellar3.analysis.differential.dge import DifferentialExpression
from src.cellar3.analysis.differential.subset import ExpressionMatrix
from src.cellar3.datasets.dataset import Dataset
from src.cellar3.export.views import ExportView
from src.cellar3.receipt import Receipt

logger = logging.getLogger('cellar.analysis.differential.exporter')


class DiffExpressionExportView(ExportView):
    analyzer = DifferentialExpression()

    def get_df(self, dataset_id: str, params: Dict[str, any]) -> DataFrame:
        meta = Receipt.create()
        target = params.get('target')
        background = params.get('background')
        dataset = Dataset(dataset_id)
        adata = ExpressionMatrix.extract(dataset=dataset, target=target, background=background, meta=meta)
        dge = DifferentialExpression().process(adata, target=target, background=background, meta=meta)
        dge = dge.fillna('')
        return dge

    def get_file_name(self, params: Dict[str, any]) -> str:
        return 'diff_gene_expression'
