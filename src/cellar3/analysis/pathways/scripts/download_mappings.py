"""
Usage:
    download_mappings.py
"""
import logging

from src.cellar3.analysis.pathways.parsers.pathways import Pathways

if __name__ == '__main__':
    logger = logging.getLogger('cellar')
    logger.addHandler(logging.StreamHandler())
    logger.setLevel(logging.INFO)
    Pathways.download_and_save()
