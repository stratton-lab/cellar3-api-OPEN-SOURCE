from typing import Dict, List

from pandas import DataFrame


class SignificantGenes:

    def __init__(self, n=5, max_pval=0.005, min_log2c=5):
        """

        :param n:
        :param max_pval: Maximum p-value for a gene to be considered.
        :param min_log2c: Minimum log fold for a gene to be considered.
        """
        self.n = n
        self.min_log2c = min_log2c
        self.max_pval = max_pval

    def extract(self, dge: DataFrame) -> List[Dict]:
        """
        Returns top genes that ar the most significant, and satisfy to minimum thresholds of pval and log2fc.
        @fixme: Make sure that we dont return genes that are not displayed (filtered) on the volcano plot.
        :param dge: Differential Gene Expression data
        :return:
        """
        significant_genes = dge[(dge['pvals'] < self.max_pval) & (abs(dge['log2_fold_change']) > self.min_log2c)]
        significant_genes = significant_genes.sort_values(by='pvals').head(self.n)
        data = significant_genes[['gene_names', 'log2_fold_change', 'minus_log10_pvals']].to_dict(orient='records')
        remap = {'gene_names': 'gene', 'log2_fold_change': 'log2fc', 'minus_log10_pvals': 'minus_log10_pval'}
        return [{remap.get(k, k): v for k, v in d.items()} for d in data]
