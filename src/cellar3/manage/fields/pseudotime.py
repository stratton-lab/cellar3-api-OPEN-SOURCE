from src.cellar3.manage.fields.smart_json import SmartJSONField
from src.cellar3.manage.widgets.list_dicts import ListOfDicts

PSEUDOTIME_FIELD = SmartJSONField(
    widget=ListOfDicts(fields=['name', 'embedding', 'directions', 'cellTypes', 'origin']),
    json_schema={
        "type": "array",
        "items": {
            "type": "object",
            "properties": {
                "name": {"type": "string"},
                "embedding": {"type": "string"},
                "cellTypes": {"type": "array", "items": {"type": "string"}},
                "directions": {
                    "type": "array",
                    "items": {
                        "type": "object",
                        "properties": {
                            "name": {"type": "string"},
                            "key": {"type": "string"}
                        },
                        "required": ["name", "key"]
                    }
                },
                "origin": {"type": "string"}
            },
            "required": ["name", "embedding", "directions"]
        }
    },
    required=False
)
