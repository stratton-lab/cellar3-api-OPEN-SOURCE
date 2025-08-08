import json
import os
import sys

import django
from django.conf import settings

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'myproject.settings')
django.setup()
print('Django SETUP')

from src.cellar3.models import DatasetMeta

if __name__ == '__main__':
    model_fields = {field.name for field in DatasetMeta._meta.get_fields()}

    with open(settings.DATASETS_META_FILE, 'r') as file:
        json_data = json.load(file)
        for item in json_data:
            filtered_item = {k: v for k, v in item.items() if k in model_fields}
            extra_fields = set(item.keys()) - model_fields
            extra_fields.discard('$schema')
            if extra_fields:
                print(f"Warning: Ignoring unexpected fields {extra_fields} in item: {item}")
            DatasetMeta.objects.create(**filtered_item)
    print('Successfully loaded data')
