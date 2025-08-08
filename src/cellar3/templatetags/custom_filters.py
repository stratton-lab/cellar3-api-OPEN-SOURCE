from django import template

register = template.Library()


@register.filter
def get_item(dictionary, key):
    if not isinstance(dictionary, dict):
        return ''
    value = dictionary.get(key)
    return value if value is not None and value != "" else ''
