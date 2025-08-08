import json
import logging
from typing import List, Any, Dict, Union

from django import forms
from django.core.exceptions import ValidationError

logger = logging.getLogger('cellar.widgets.list_dict')


class ListOfDicts(forms.Widget):
    """
    Maps custom keys to custom values.
    Only supports strings for keys and values.
    """
    template_name = 'widgets/json_table_widget.html'

    def __init__(self, *args, **kwargs):
        self.fields = kwargs.pop('fields', [])
        self.json_schema = kwargs.pop('json_schema', None)
        self.debug_watch_fields = {}
        super().__init__(*args, **kwargs)

    def _load_json(self, value: Union[str, List[Dict]]) -> List[Dict]:
        """
        Loads the json sting stored in a DB cell as a python list of dicts.
        :param value:
        :return:
        """
        if isinstance(value, str):
            try:
                return json.loads(value)
            except json.JSONDecodeError:
                return []
        # If already List
        return value

    def format_value(self, value: str) -> List[Dict]:
        """
        Used to decide how to display a DB cell
        :param value:
        :return:
        """
        rows = self._load_json(value)

        if not rows:  # When creating a new Dataset, the value is None
            return []

        if not isinstance(rows, list):
            raise ValidationError(f'Incorrect JSON in ListOfDicts field. Expected List, got {rows}')
        self._convert_database2display(rows)
        return rows

    def value_from_datadict(self, data, files, name):
        """
        Converting into the object to save in DB.
        :param data:
        :param files:
        :param name:
        :return:
        """
        list_dicts = []

        for row in zip(*self._get_field_values(data, name)):
            # If field declared in field_types, we convert the value to right type
            row_dict = {field: self._convert_display2database(value, field) for field, value in zip(self.fields, row)}
            # Avoiding saving rows if Name (first cell) empty
            if row_dict.get(self.fields[0], None) and row_dict:
                list_dicts.append(row_dict)

        value = json.dumps(list_dicts)
        if name in self.debug_watch_fields:
            logger.info(f'[ListOfDicts] exporting field {name} TYPE {type(value)} as {value}')

        return value

    def _get_expected_field_type(self, field: str) -> str:
        """
        Returns:
        - json
        - string
        - number
        :param field:
        :return:
        """
        if self.json_schema:
            expected_type = self.json_schema["items"]["properties"][field]["type"]
            if expected_type in {"array", "object"}:
                return 'json'
            return expected_type
        return 'string'  # If no schema provided, all values treated as strings

    def _convert_display2database(self, value: Any, field: str):
        if self.json_schema:
            try:
                if value and self._get_expected_field_type(field) == 'json':
                    return json.loads(value)
            except Exception as e:
                logger.error(f'Error converting value {value} from field {field} to JSON: {e}')
        return value

    def _convert_database2display(self, rows: List[Dict]):
        if self.json_schema:
            for row in rows:
                for field, value in row.items():
                    if self._get_expected_field_type(field) == 'json':
                        row[field] = json.dumps(value)

    def _get_field_values(self, data, name: str) -> List:
        field_values = []
        for field in self.fields:
            field_data = data.getlist(f'{name}_{field}')
            field_values.append(field_data)
        return field_values

    def get_context(self, name, value, attrs):
        """
        Used for displaying the field in admin interface.
        :param name:
        :param value:
        :param attrs:
        :return:
        """
        context = super().get_context(name, self.format_value(value), attrs)
        context['widget']['fields'] = self.fields  # List of dict keys from main dict
        context['widget']['value'] = self.format_value(value)  # List of dicts
        return context
