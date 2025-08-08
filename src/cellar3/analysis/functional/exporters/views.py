import logging
from typing import Dict

import numpy as np
from pandas import DataFrame

from src.cellar3.analysis.differential.cached import get_cached_dge
from src.cellar3.analysis.functional.enrich import Enrichir
from src.cellar3.export.views import ExportView
from src.cellar3.receipt import Receipt

logger = logging.getLogger('cellar.analysis.functional.exporter')


class FunctionalExportView(ExportView):
    analyzer = Enrichir()

    def get_df(self, dataset_id: str, params: Dict[str, any]) -> DataFrame:
        meta = Receipt.create()
        target = params.get('target')
        background = params.get('background')
        gene_set = params.get('gene_set')
        diff_p_val_thr = 0.05  # P-value Threshold
        fc_thr = 0.25  # Fold Change Threshold
        dataset, meta, dge = get_cached_dge(dataset_id, target, background, diff_p_val_thr, fc_thr, meta)
        df = self.analyzer.process(gene_set=gene_set, dge=dge, organism=dataset.meta['species'], meta=meta)
        df['minus_log10_pvals'] = -np.log10(df['Adjusted P-value'])
        term2genes = {term: Enrichir.get_dge_genes_in_pathway(dge, df, term) for term in df['Term'].unique()}

        exported_df = df[['Gene_set', 'Term', 'minus_log10_pvals', 'Genes', 'Overlap']].copy()
        self._add_dge_column(exported_df, term2genes, 'upregulated')
        self._add_dge_column(exported_df, term2genes, 'downregulated')

        return exported_df

    @classmethod
    def _add_dge_column(cls, df: DataFrame, term2genes: Dict, column: str):
        df[column] = df['Term'].apply(lambda term: ';'.join(term2genes[term][column]) if term in term2genes else [])

    def get_file_name(self, params: Dict[str, any]) -> str:
        return 'functional_analysis'
