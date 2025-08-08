from src.cellar3.manage.fields.smart_json import SmartJSONField
from src.cellar3.manage.widgets.list_dicts import ListOfDicts

LINKS_FIELD = SmartJSONField(
    widget=ListOfDicts(fields=['name', 'url']),
    json_schema={
        "type": "array",
        "items": {
            "type": "object",
            "properties": {
                "name": {
                    "type": "string",
                    "minLength": 1
                },
                "url": {
                    "type": "string",
                    "minLength": 1,
                    "format": "uri"
                }
            },
            "required": ["name", "url"]
        }
    },
    required=False
)
