import logging
from time import sleep
from typing import Dict, List
from xml.etree import ElementPath

import requests

from src.cellar3.analysis.pathways.exceptions import GeneException
from src.cellar3.analysis.pathways.parsers.entry import EntryParser

logger = logging.getLogger('cellar.analysis.pathways.gene')


class Gene(EntryParser):
    e_type = 'gene'

    @classmethod
    def match(cls, n: ElementPath):
        return n.get('type') == cls.e_type and n.find('graphics').get('type') == 'rectangle'

    @classmethod
    def parse_gene_symbols(cls, text: str) -> Dict[str, List[str]]:
        """
        Note: Some genes do not have symbols, ex: https://www.kegg.jp/dbget-bin/www_bget?hsa:107987478
        :param text:
        :return:
        """
        gene2symbols = {}
        docs = text.strip().split('///\n')
        for doc in docs:
            if doc:
                code = None
                organism = None
                symbols = []  # default if gene has no symbols
                lines = doc.split('\n')
                for line in lines:
                    if line.startswith('ENTRY'):
                        code = line.replace('ENTRY', '').strip().split(' ')[0]
                    if line.startswith('ORGANISM'):
                        organism = line.replace('ORGANISM', '').strip().split(' ')[0]
                    if line.startswith('SYMBOL'):
                        symbols = line.replace('SYMBOL', '').strip().split(', ')
                if code and organism:
                    gene_id = f'{organism}:{code}'
                    gene2symbols[gene_id] = symbols
                else:
                    raise GeneException(f'Missing ENTRY or ORGANISM field values in KEGG doc: {doc}')

        return gene2symbols

    @staticmethod
    def chunk_list(lst, chunk_size=10):
        return [lst[i:i + chunk_size] for i in range(0, len(lst), chunk_size)]

    @classmethod
    def qc(cls, gene_ids: List[str], gene2symbols: Dict[str, List[str]]):
        if len(gene2symbols.keys()) != len(gene_ids):
            missing_genes = set(gene_ids) - set(gene2symbols.keys())
            raise GeneException(
                f'Generated gene2symbols map for {len(gene2symbols.keys())}/{len(gene_ids)}. Missing: {missing_genes}')

    @classmethod
    def download_gene2symbols(cls, gene_ids: List[str]) -> Dict[str, List[str]]:
        url = f'https://rest.kegg.jp/get/{"+".join(gene_ids)}'
        res = requests.get(url)
        res.raise_for_status()
        gene2symbols = cls.parse_gene_symbols(res.text)
        cls.qc(gene_ids, gene2symbols)
        return gene2symbols

    @classmethod
    def download_gene2symbols_chunks(cls, all_gene_ids: List[str]) -> Dict[str, List[str]]:
        gene2symbols = {}
        for gene_ids in cls.chunk_list(all_gene_ids):
            gene2symbols.update(cls.download_gene2symbols(gene_ids))
            sleep(0.33)
        return gene2symbols

    @classmethod
    def assign_symbols(cls, gene: Dict, gene2symbols: Dict[str, List[str]]):
        for gene_id in gene.get('gene_ids', []):
            symbols = gene2symbols.get(gene_id)
            if symbols is None:
                raise GeneException(f'Missing symbols for gene {gene_id}')
            gene['symbols'].extend(symbols)

    @classmethod
    def collect_gene_ids(cls, genes: List[Dict]) -> List[str]:
        gene_ids = set()
        for gene in genes:
            gene_ids.update(gene['gene_ids'])
        return list(gene_ids)

    @classmethod
    def collect_and_assign_symbols(cls, genes: List[Dict]):
        gene_ids = Gene.collect_gene_ids(genes)
        gene2symbols = Gene.download_gene2symbols_chunks(gene_ids)
        for gene in genes:
            cls.assign_symbols(gene, gene2symbols)

    @classmethod
    def parse(cls, n: ElementPath) -> Dict:
        gene_id = n.get('id')
        try:
            x = cls.graph_int(n, 'x')
            y = cls.graph_int(n, 'y')
            width = cls.graph_int(n, 'width')
            height = cls.graph_int(n, 'height')
            return {
                'id': gene_id,
                'x': x - (width / 2),
                'y': y - (height / 2),
                'width': width,
                'height': height,
                'gene_ids': n.get('name').split(' '),
                'symbols': []
            }
        except Exception as e:
            raise GeneException(f'Could not parse gene {gene_id} : {e}')
