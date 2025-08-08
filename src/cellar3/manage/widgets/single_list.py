import json
import logging
from typing import Tuple, List, Dict, Union

from django import forms

from src.cellar3.datasets.categories import CATEGORIES_READER

logger = logging.getLogger('cellar.widgets.single_dict')


class SingleValueJSONListWidget(forms.Select):

    def __init__(self, *args, **kwargs):
        choices = kwargs.pop('choices', [])
        super().__init__(*args, **kwargs)
        self.choices = choices
        self.debug = False

    def value_from_datadict(self, data, files, name):
        raw_value = data.get(name)
        value = [raw_value] if raw_value else []
        if self.debug:
            logger.info(f'[SingleValueJSONListWidget] exporting field {name} as {value}')
        return json.dumps(value)

    def format_value(self, value):
        if isinstance(value, str):
            try:
                value = json.loads(value)
            except json.JSONDecodeError:
                pass
        if isinstance(value, list) and len(value) == 1:
            return value[0]
        return super().format_value(value)


class CategoriesListWidget(SingleValueJSONListWidget):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # @todo Load before display instead of at server start
        self.choices = self.load_choices()

    @classmethod
    def load_choices(cls) -> List[Tuple[str, str]]:
        categories_root = CATEGORIES_READER.get_categories()
        choices = []
        cls.categories2choices(categories_root, None, choices)
        return choices

    @classmethod
    def categories2choices(cls, category: Dict, name_prefix: Union[str, None], choices: List[Tuple[str, str]]):
        category_id = category.get('id', '-')
        category_name = category.get('name', '')
        category_full_name = ' - '.join([name_prefix, category_name]) if name_prefix else category_name
        choices.append((category_id, category_full_name))

        for subcategory in category.get('categories', []):
            cls.categories2choices(subcategory, category_full_name, choices)
