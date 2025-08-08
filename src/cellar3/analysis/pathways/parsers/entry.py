from typing import Dict, List
from xml.etree import ElementPath

from src.cellar3.analysis.pathways.exceptions import ParsingException


class EntryParser:
    e_type = None

    @classmethod
    def match(cls, n: ElementPath):
        return n.get('type') == cls.e_type

    @classmethod
    def parse(cls, n: ElementPath) -> Dict:
        pass

    @classmethod
    def parse_all(cls, root: ElementPath) -> List[Dict]:
        """
        :param root:
        :return:
        """
        elements = []
        for n in root.findall('entry'):
            if cls.match(n):
                elements.append(cls.parse(n))
        return elements

    @classmethod
    def graph_int(cls, n: ElementPath, field: str) -> int:
        """
        Helper function. Returns a XML attribute cast to integer.
        :param n:
        :param field:
        :return:
        """
        try:
            return int(n.find('graphics').get(field))
        except TypeError as e:
            raise ParsingException(f'Could not cast to integer the field {field}')
