import logging
from typing import Dict, List

from src.cellar3.analysis.pathways.parsers.pathway import Pathway
from src.cellar3.analysis.pathways.species2code import SPECIES2CODE
from src.cellar3.analysis.pathways.term2id import TERM2ID
from src.cellar3.exceptions import UserException

logger = logging.getLogger('cellar.analysis.pathways.kegg')


class KEGG:

    def __init__(self):
        self.term2pathway = TERM2ID
        self.pathways_ids = set(self.term2pathway.values())
        self.species2code = SPECIES2CODE

    """
    Term can be a pathway name or pathway id
    """
    def term2id(self, term: str) -> str:
        if term in self.pathways_ids:
            return term
        return self.term2pathway.get(term)

    def get_pathway(self, species: str, pathway: str, upregulated: List[str], downregulated: List[str]) -> Dict:
        species_code = self.species2code.get(species)
        if not species_code:
            raise UserException(f'Species not available: {species}')
        pathway_map_id = self.term2id(pathway)
        if not pathway_map_id:
            raise UserException(f'No pathway maps for term {pathway}')

        pathway = Pathway.load(species_code, pathway_map_id)
        genes = pathway.get('genes', [])
        self.mark_dge_genes(genes, upregulated, downregulated)

        return pathway

    @classmethod
    def mark_dge_genes(cls, genes: List[Dict], upregulated: List[str], downregulated: List[str]):
        """
        Sets status of down or upregulated genes.
        @todo log exception if down/up gene not found in pathway
        :param genes:
        :param upregulated:
        :param downregulated:
        :return:
        """
        status2genes = {'upregulated': set(upregulated), 'downregulated': set(downregulated)}
        for gene in genes:
            for symbol in gene.get('symbols', []):
                for status, status_genes in status2genes.items():
                    if symbol in status_genes:
                        gene['status'] = status
                        gene['symbols'] = [symbol]  # Displaying matched symbol as gene name
