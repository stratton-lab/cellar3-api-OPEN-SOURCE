import json
import logging

from django import forms

logger = logging.getLogger('cellar.widgets.def_dict')


class PredefinedKeysJSONWidget(forms.Widget):
    """
    Uses a mapping dict with predefined keys. Used for fields like Maintainer.
    The map has to be a simple, one level dict as {key:value}.
    Value is string.
    """
    template_name = 'widgets/predefined_keys_json_widget.html'

    def __init__(self, *args, **kwargs):
        self.predefined_keys = kwargs.pop('predefined_keys', [])
        super().__init__(*args, **kwargs)

    def format_value(self, value):
        if isinstance(value, str):
            try:
                value = json.loads(value)
            except json.JSONDecodeError:
                value = {}
        return value

    def value_from_datadict(self, data, files, name):
        """
        Only exporting values that are not null and not None
        :param data:
        :param files:
        :param name:
        :return:
        """
        dict_data = {}
        for key in self.predefined_keys:
            value = data.get(f'{name}_{key}', '')
            if value:
                dict_data[key] = value
        # logger.info(f'[PredefinedKeysJSONWidget] exporting field {name} as {dict_data}')
        return json.dumps(dict_data)

    def get_context(self, name, value, attrs):
        context = super().get_context(name, self.format_value(value), attrs)
        context['widget'].update({
            'predefined_keys': self.predefined_keys,
            'value': self.format_value(value)
        })
        return context
