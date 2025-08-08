from src.cellar3.manage.fields.smart_json import SmartJSONField
from src.cellar3.manage.widgets.dict import JSONDictWidget

INFO_FIELD = SmartJSONField(
    widget=JSONDictWidget(),
    json_schema={
        "type": "object",
        "propertyNames": {
            "minLength": 1
        },
        "additionalProperties": {
            "type": "string",
            "minLength": 1
        }
    },
    required=False
)
