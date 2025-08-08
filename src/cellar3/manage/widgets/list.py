import json
import logging

from django import forms

logger = logging.getLogger('cellar.widgets.list')


class ListWidget(forms.TextInput):
    """
    A custom widget that displays a JSON list as a comma-separated string
    and parses it back into a JSON list.
    """

    def __init__(self, *args, **kwargs):
        self.separator = kwargs.pop('separator', ',')
        super().__init__(*args, **kwargs)

    def format_value(self, value):
        """
        Format the value for display in the widget.
        """
        if isinstance(value, str):
            # Try to load the JSON string if the value is a string
            try:
                value = json.loads(value)
            except json.JSONDecodeError:
                value = []
        if isinstance(value, list):
            # Convert list to a comma-separated string
            return f'{self.separator} '.join(map(str, value))
        return super().format_value(value)

    def value_from_datadict(self, data, files, name):
        """
        Parse the input from the form back into a JSON list.
        """
        value = data.get(name)
        if value:
            return json.dumps([item.strip() for item in value.split(self.separator) if item.strip()])
        return json.dumps([])
