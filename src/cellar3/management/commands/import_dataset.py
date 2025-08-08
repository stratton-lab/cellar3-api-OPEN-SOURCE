import json

from django.core.management import BaseCommand

from src.cellar3.models import DatasetMeta


class Command(BaseCommand):
    help = 'Imports a dataset in JSON format into the DB.'

    def add_arguments(self, parser):
        parser.add_argument('file_path', type=str, help='Path to the JSON file.')

    def handle(self, *args, **kwargs):
        try:
            file_path = kwargs['file_path']
            with open(file_path, 'r') as file:
                json_data = json.load(file)
                unique_fields = {"id": json_data["id"]}
                fields_to_update = {key: value for key, value in json_data.items() if key != "id"}
                obj, created = DatasetMeta.objects.update_or_create(defaults=fields_to_update, **unique_fields)
                status = "CREATED" if created else "UPDATED"
            self.stdout.write(self.style.SUCCESS(f'Successfully imported dataset from {file_path}: {status}.'))
        except BaseException as e:
            self.stdout.write(self.style.ERROR(str(e)))
