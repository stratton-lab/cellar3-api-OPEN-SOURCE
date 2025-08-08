import os

from django.conf import settings
from django.core.management import BaseCommand

from src.cellar3.models import DatasetMeta


class Command(BaseCommand):
    help = 'Deletes a dataset (meta record in DB and h5ad file on disk.'

    def add_arguments(self, parser):
        # Optional: Define command-line arguments here
        parser.add_argument('dataset_id', type=str, help='ID of the dataset to delete.')

    def handle(self, *args, **kwargs):
        try:
            dataset_id = kwargs['dataset_id']
            dataset = DatasetMeta.objects.get(pk=dataset_id)
            self._safe_delete_dataset_file(dataset.file)
            dataset.delete()
            self.stdout.write(self.style.SUCCESS(f'Dataset deleted for id {dataset_id}'))
        except Exception as e:
            self.stdout.write(self.style.ERROR(e))

    @classmethod
    def _safe_delete_dataset_file(cls, file_name: str):
        file_path = os.path.join(settings.DATASETS_FOLDER, file_name)
        if file_path and os.path.exists(file_path):
            os.remove(file_path)
