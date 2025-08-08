import base64
import json
import logging
import os.path
import struct
import xml.etree.ElementTree as ET
from time import sleep
from typing import Dict, List, Set

import requests
from django.conf import settings
from tqdm import tqdm

from src.cellar3.analysis.pathways.exceptions import PathwayUnavailableException
from src.cellar3.analysis.pathways.parsers.gene import Gene
from src.cellar3.exceptions import UserException

logger = logging.getLogger('cellar.analysis.pathways.pathway')

IGNORE_PATHWAYS = {
    '01100'  # Metabolic Pathways, huge representation of all pathways.
}



class Pathway:

    @classmethod
    def load(cls, species_code: str, pathway_id: str) -> Dict:
        """
        Loads a previously downloaded and saved pathway.
        :param species_code:
        :param pathway_id:
        :return:
        """
        if cls.has_saved_pathway(species_code, pathway_id):
            with open(cls.get_pathway_path(species_code, pathway_id), 'r') as file:
                return json.load(file)
        raise UserException(f'Pathway {pathway_id} for species {species_code} is not yet available.')

    @classmethod
    def parse(cls, pathway_xml: str) -> Dict:
        """
        Returns a dict representation of a pathway with a list of gene nodes.
        :param pathway_xml:
        :return:
        """
        root = ET.fromstring(pathway_xml)
        pathway_id = root.get('number')
        species = root.get('org')
        try:
            genes = Gene.parse_all(root)
            Gene.collect_and_assign_symbols(genes)
            img_base64, width, height = cls.download_image(root.get('image'))
            return {
                'id': pathway_id,
                'species': species,
                'name': root.get('title'),
                'width': width,
                'height': height,
                'genes': genes,
                'base64': img_base64
            }
        except Exception as e:
            raise UserException(f'Could not parse pathway {pathway_id} for species {species} : {e}')

    @classmethod
    def download_image(cls, url: str):
        res = requests.get(url)
        res.raise_for_status()
        if res.content[:8] != b'\x89PNG\r\n\x1a\n':
            raise UserException("Not a valid PNG file")
        img_base64 = base64.b64encode(res.content).decode('utf-8')
        width, height = struct.unpack('>II', res.content[16:24])
        return img_base64, width, height

    @classmethod
    def download(cls, species_code: str, pathway_id: str) -> Dict:
        """
        Downloads the basic information about a pathway and a list of genes with info on how to display them.
        Ex: https://rest.kegg.jp/get/eco00020/kgml
        :param species_code:
        :param pathway_id:
        :return:
        """
        full_pathway_id = f'{species_code}{pathway_id}'
        url = f'http://rest.kegg.jp/get/{full_pathway_id}/kgml'
        res = requests.get(url)
        if res.status_code == 404:
            raise PathwayUnavailableException(f'Pathway {pathway_id} not available for species {species_code}.')
        res.raise_for_status()
        pathway_xml = res.text.strip()
        return Pathway.parse(pathway_xml)

    @classmethod
    def has_saved_pathway(cls, species_code: str, pathway_id: str) -> bool:
        return os.path.exists(cls.get_pathway_path(species_code, pathway_id))

    @classmethod
    def get_pathway_path(cls, species_code: str, pathway_id: str) -> str:
        return f'{settings.PATHWAYS_FOLDER}/{species_code}/{species_code}{pathway_id}.json'

    @classmethod
    def get_ignored_fle_path(cls, species_code: str):
        return f'{settings.PATHWAYS_FOLDER}/{species_code}/ignored.tsv'

    @classmethod
    def add_ignored_pathway(cls, species_code: str, pathway_id: str):
        with open(cls.get_ignored_fle_path(species_code), 'a') as file:
            file.write(f'{pathway_id}\n')

    @classmethod
    def load_ignored_pathways(cls, species_code: str) -> Set[str]:
        ignored = IGNORE_PATHWAYS
        ignored_path = cls.get_ignored_fle_path(species_code)
        if os.path.exists(ignored_path):
            with open(ignored_path, 'r') as file:
                for line in file.readlines():
                    ignored.add(line.strip())
        return ignored

    @classmethod
    def download_and_save(cls, species_code: str, pathway_ids: List[str]):
        """
        :param species_code:
        :param pathway_ids:
        :return:
        """
        ignored = cls.load_ignored_pathways(species_code)
        for pathway_id in tqdm(pathway_ids):
            if pathway_id not in ignored:
                path = cls.get_pathway_path(species_code, pathway_id)
                if not cls.has_saved_pathway(species_code, pathway_id):
                    try:
                        pathway = cls.download(species_code, pathway_id)
                        with open(path, 'w') as file:
                            json.dump(pathway, file)
                    except PathwayUnavailableException as e:
                        logger.info(f'Skipping : {e}')
                        cls.add_ignored_pathway(species_code, pathway_id)
                    sleep(0.3)
                else:
                    logger.info(f'Skipping previously downloaded pathway {species_code} {pathway_id}')
            else:
                logger.info(f'Ignoring pathway {pathway_id} as it is on the ignore list.')
