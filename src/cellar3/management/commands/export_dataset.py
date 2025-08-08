import json

from django.core.management import BaseCommand

from src.cellar3.models import DatasetMeta
from src.cellar3.serializers import DatasetInfoSerializer


class Command(BaseCommand):
    help = 'Exports a dataset row from the DB into a JSON file.'

    def add_arguments(self, parser):
        parser.add_argument('dataset_id', type=str, help='ID of the dataset to export.')

    def handle(self, *args, **kwargs):
        try:
            dataset_id = kwargs['dataset_id']
            dataset_meta = DatasetMeta.objects.get(pk=dataset_id)
            serializer = DatasetInfoSerializer(dataset_meta)
            row = serializer.data
            self.stdout.write(json.dumps(row))
            # self.stdout.write(self.style.SUCCESS(f'Successfully exported dataset {dataset_id}.'))
        except BaseException as e:
            self.stdout.write(self.style.ERROR(str(e)))
