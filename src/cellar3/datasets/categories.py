import logging
import os
from typing import Any, Dict

import yaml
from django.conf import settings

logger = logging.getLogger('cellar.categories')


class AbstractCategoriesReader:
    def get_categories(self) -> Dict[str, Any]:
        pass


class YamlCategoriesReader(AbstractCategoriesReader):

    def get_categories(self) -> Dict[str, Any]:
        if os.path.exists(settings.CATEGORIES_FILE):
            return yaml.safe_load(open(settings.CATEGORIES_FILE, 'r'))
        else:
            logger.warning(f'Could not load categories from {settings.CATEGORIES_FILE}. Returning empty dictionary.')
            return {}


CATEGORIES_READER = YamlCategoriesReader()
