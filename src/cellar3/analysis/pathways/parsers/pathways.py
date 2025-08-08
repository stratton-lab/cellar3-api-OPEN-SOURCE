import json
import logging
from typing import Dict

import requests

from src import settings

logger = logging.getLogger('cellar.analysis.pathways.mapping')


class Pathways:

    @classmethod
    def get_pathways_mappings_path(cls) -> str:
        return f'{settings.PATHWAYS_FOLDER}/term2id.json'

    @staticmethod
    def download_pathways_mapping() -> Dict:
        """
        Maps terms (pathway names) to patheay map ids (numbers)
        The downloaded mapping is not species specific. i.e. some pathways may not appear in some species.
        From https://rest.kegg.jp/list/pathway
        (species specific:         )
        :return:
        """
        url = 'https://rest.kegg.jp/list/pathway'
        res = requests.get(url)
        res.raise_for_status()
        lines = res.text.strip().split('\n')
        mapping = {}
        for line in lines:
            if line:
                parts = line.split('\t')
                if len(parts) != 2:
                    raise Exception(f'Invalid line in pathways list file: {line}')
                code, term = parts
                mapping[term] = code.replace('map', '')
        return mapping

    @classmethod
    def download_and_save(cls):
        path = cls.get_pathways_mappings_path()
        mapping = cls.download_pathways_mapping()
        with open(path, 'w') as file:
            json.dump(mapping, file)
