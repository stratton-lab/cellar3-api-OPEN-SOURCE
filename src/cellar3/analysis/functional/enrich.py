import abc
import logging
import warnings
from typing import Dict, List, Set

import gseapy as gp
from pandas import DataFrame

from src.cellar3.analysis.differential.dge import DifferentialExpression
from src.cellar3.exceptions import UserException
from src.cellar3.receipt import Receipt
from src.cellar3.tools import timeit

logger = logging.getLogger('cellar.analysis.expression.enrichment')

FIELD_TERM = 'Term'
FIELD_MINUS_LOG_10P = 'minus_log10_pvals'


class EnrichmentAnalysis:
    @abc.abstractmethod
    def process(self, gene_set: str, dge: DataFrame, organism: str, meta: Dict[str, any]) -> DataFrame:
        pass


class Enrichir(EnrichmentAnalysis):
    """
    Cellar1 code for enrichment: https://github.com/euxhenh/cellar/blob/a1b9e551ed1ba29a7808df9dc9cc4359aae27008/controller/cellar/core/_de.py#L97
    GSEAPY library: https://gseapy.readthedocs.io/en/latest/introduction.html#gseapy-gene-set-enrichment-analysis-in-python

    """

    def __init__(self, cutoff=0.05):
        self.cutoff = cutoff

    @timeit
    def process(self, gene_set: str, dge: DataFrame, organism: str, meta: Dict[str, any]) -> DataFrame:
        """

        :param dge:
        :param gene_set:
        :param organism:
        :param meta:
        :return:
        """
        upregulated = DifferentialExpression.get_genes(dge, DifferentialExpression.CATEGORY_UP_REGULATED)
        downregulated = DifferentialExpression.get_genes(dge, DifferentialExpression.CATEGORY_DOWN_REGULATED)
        # Gene symbols are all “UPPERCASE” in the Enrichr Libraries.
        target = [gene.upper() for gene in upregulated + downregulated]
        if not target:
            raise UserException('No genes are significantly up or down regulated.')

        logger.info(f'Enriching {len(target)} genes with gene_set {gene_set}')
        warnings.filterwarnings('ignore', category=FutureWarning, module='pandas')
        Receipt.add_enrichir_step(meta, params={'organism': organism, 'cutoff': self.cutoff, 'gene_sets': gene_set})
        # logger.info(f'Available gene sets: {gp.get_library_name()}')
        enr = gp.enrichr(
            organism=organism,
            gene_list=target,
            gene_sets=gene_set,
            outdir=None,
            no_plot=True,
            cutoff=self.cutoff
        )
        return enr.res2d

    @staticmethod
    def filter_genes(genes: List[str], matching: Set[str]) -> List[str]:
        """
        Only keeps genes from genes if their uppercased version is also in matching.
        :param genes:
        :param matching:
        :return:
        """
        return [gene for gene in genes if gene.upper() in matching]

    @staticmethod
    def get_pathway_genes(enrichment: DataFrame, pathway: str) -> Set[str]:
        """
        Returns all genes in a given pathway, corresponding to an enriched term.
        :param enrichment:
        :param pathway:
        :return:
        """
        filtered_df = enrichment[enrichment['Term'] == pathway]
        return set(filtered_df["Genes"].iloc[0].split(';'))

    @staticmethod
    def get_pathway_total(enrichment: DataFrame, pathway: str) -> int:
        """
        Returns the total number of genes in that pathway
        :param enrichment:
        :param pathway:
        :return:
        """
        try:
            filtered_df = enrichment[enrichment['Term'] == pathway]
            return int(filtered_df["Overlap"].iloc[0].split('/')[1])
        except Exception as e:
            logger.error(f'Could not parse term overlap: {e}')
            return -1

    @staticmethod
    def get_dge_genes_in_pathway(dge: DataFrame, enrichment: DataFrame, pathway: str) -> Dict[str, List[str]]:
        """

        :param dge:
        :param enrichment:
        :param pathway:
        :return:
        """
        matched_genes = Enrichir.get_pathway_genes(enrichment, pathway)
        upregulated = Enrichir.filter_genes(DifferentialExpression.get_up_genes(dge), matched_genes)
        downregulated = Enrichir.filter_genes(DifferentialExpression.get_down_genes(dge), matched_genes)
        total = Enrichir.get_pathway_total(enrichment, pathway)
        return {'upregulated': upregulated, 'downregulated': downregulated, 'total': total}

    @staticmethod
    def get_term_genes(dge: DataFrame, enrichment: DataFrame, top_terms=10) -> Dict[str, Dict]:
        terms = enrichment.sort_values(by=FIELD_MINUS_LOG_10P, ascending=False).head(top_terms)[FIELD_TERM].unique()
        return {term: Enrichir.get_dge_genes_in_pathway(dge, enrichment, term) for term in terms}

    @staticmethod
    def get_library(pathway: str, species: str):
        return gp.get_library(name=pathway, organism=species)


class Prerank(EnrichmentAnalysis):
    def process(self, gene_set: str, dge: DataFrame, organism: str, meta: Dict[str, any]) -> DataFrame:
        """
        @todo Implement to use for Network vizu
        https://gseapy.readthedocs.io/en/latest/gseapy_example.html#Prerank-example
        :return:
        """
        dge['rank_metric'] = dge['log2_fold_change'] * dge['minus_log10_pvals']
        dge['gene_names'] = dge['gene_names'].str.upper()
        df_sorted = dge.sort_values(by='rank_metric', ascending=False)
        ranked_list = df_sorted[['gene_names', 'rank_metric']]
        logger.info(f'Enriching {ranked_list} genes [PREPRANK] with gene_set {gene_set}')
        # FIXME: KeyError gene_names ( doesn't always happens)
        enr = gp.prerank(
            rnk=ranked_list,
            gene_sets=gene_set,
            outdir=None,
            no_plot=True,
            min_size=5,
            max_size=1000,
            seed=6,
            verbose=True
        )
        return enr.res2d
