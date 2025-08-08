from cellar3.manage.fields.smart_json import SmartJSONField
from cellar3.manage.widgets.list_dicts import ListOfDicts

# TODO Potentially empty list of dicts
GROUPS_FIELD = SmartJSONField(
    widget=ListOfDicts(fields=['name', 'key']),
    json_schema={
        "type": "array",
        "items": {
            "type": "object",
            "properties": {
                "name": {
                    "type": "string",
                    "minLength": 1
                },
                "key": {
                    "type": "string",
                    "minLength": 1
                }
            },
            "required": ["name", "key"],
            "additionalProperties": False
        },
        "minItems": 0
    },
    required=True
)
