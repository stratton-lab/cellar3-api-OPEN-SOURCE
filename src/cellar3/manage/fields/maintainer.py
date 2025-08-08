from src.cellar3.manage.fields.smart_json import SmartJSONField
from src.cellar3.manage.widgets.def_dict import PredefinedKeysJSONWidget

MAINTAINER_FIELD = SmartJSONField(
    widget=PredefinedKeysJSONWidget(predefined_keys=['name', 'email', 'affiliation']),
    json_schema={
        "type": "object",
        "properties": {
            "name": {
                "type": "string",
                "minLength": 1
            },
            "email": {
                "type": "string",
                "minLength": 1,
                "format": "email"
            },
            "affiliation": {
                "type": "string",
                "minLength": 1
            }
        },
        "required": ["name", "email", "affiliation"]
    },
    required=False
)
