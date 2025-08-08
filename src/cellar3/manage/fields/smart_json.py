import logging

from django.core.exceptions import ValidationError
from django.forms import JSONField
from jsonschema import ValidationError as JSONSchemaValidationError
from jsonschema.validators import validate as jsonschema_validate, Draft202012Validator

logger = logging.getLogger('cellar.fields.smart_json')


class SmartJSONField(JSONField):
    # Having [] in default empty_values list prevented an empty list [] to be saved into the DB.
    empty_values = [None, ""]

    def __init__(self, json_schema=None, widget=None, *args, **kwargs):
        self.json_schema = json_schema
        widget.json_schema = json_schema
        super().__init__(widget=widget, *args, **kwargs)

    def clean(self, value):
        # First, call the parent's clean method to handle basic validation
        value = super().clean(value)

        # Now, perform custom validation using the JSON schema
        if self.json_schema:
            try:
                jsonschema_validate(instance=value, schema=self.json_schema,
                                    format_checker=Draft202012Validator.FORMAT_CHECKER)
            except JSONSchemaValidationError as e:
                raise ValidationError(f"Invalid input: {e.message}")

        # Making sure we still have JSON at the end
        if not isinstance(value, (dict, list)):
            raise ValidationError(f"Incorrect data type, expected dict or list, but got: {value.__class__.__name__}")

        return value
