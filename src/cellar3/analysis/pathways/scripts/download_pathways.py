"""
Usage:
    download_pathways.py
"""
import logging

from src.cellar3.analysis.pathways.parsers.pathway import Pathway
from src.cellar3.analysis.pathways.parsers.pathways import Pathways
from src.cellar3.analysis.pathways.parsers.species import Species

if __name__ == '__main__':
    logger = logging.getLogger('cellar')
    logger.addHandler(logging.StreamHandler())
    logger.setLevel(logging.INFO)
    for species_code in Species.load().values():
        map_ids = list(Pathways.load().values())
        Pathway.download_and_save(species_code, map_ids)
