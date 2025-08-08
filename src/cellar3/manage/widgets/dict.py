import json
import logging

from django import forms

logger = logging.getLogger('cellar.widgets.dict')


class JSONDictWidget(forms.Widget):
    """
    Displays the JSON as text
    """
    template_name = 'widgets/json_dict_widget.html'

    def __init__(self, *args, **kwargs):
        self.debug_watch_fields = {}
        super().__init__(*args, **kwargs)

    def format_value(self, value):
        if isinstance(value, str):
            try:
                value = json.loads(value)
            except json.JSONDecodeError:
                value = {}
        return value

    def value_from_datadict(self, data, files, name):
        keys = data.getlist(f'{name}_keys[]')
        values = data.getlist(f'{name}_values[]')
        dict_data = {}
        for key, value in zip(keys, values):
            if key.strip():  # only add non-empty keys
                dict_data[key] = value

        if name in self.debug_watch_fields:
            logger.info(f'[JSONDictWidget] exporting field {name} as {dict_data}')

        return json.dumps(dict_data)

    def get_context(self, name, value, attrs):
        context = super().get_context(name, self.format_value(value), attrs)
        context['widget'].update({
            'value': self.format_value(value)
        })
        return context
