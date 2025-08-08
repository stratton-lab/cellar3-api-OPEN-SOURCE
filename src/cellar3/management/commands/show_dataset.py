import anndata
from django.conf import settings
from django.core.management import BaseCommand

from src.cellar3.models import DatasetMeta

GREEN = "\033[32m"
BLUE = "\033[36m"
RESET = "\033[0m"

class Command(BaseCommand):
    help = 'Displays all meta fields about a dataset.'

    def add_arguments(self, parser):
        parser.add_argument('dataset_id', type=str, help='ID of the dataset to display')

    def handle(self, *args, **kwargs):
        try:
            dataset_id = kwargs['dataset_id']
            dataset_meta = DatasetMeta.objects.get(pk=dataset_id)
            self.stdout.write(self.style.SUCCESS(f'Displaying fields for dataset ID: {dataset_id}\n'))

            for field in dataset_meta._meta.fields:
                field_name = field.verbose_name if hasattr(field, 'verbose_name') else field.name
                field_value = getattr(dataset_meta, field.name, None)

                # Display field name in bold and blue
                self.stdout.write(f'{BLUE}{field_name.capitalize()}{RESET}: {field_value}')

            # Also display embeddings and columns
            self.stdout.write(self.style.SUCCESS('DATASET OBJECT'))
            adata = anndata.read_h5ad(f'{settings.DATASETS_FOLDER}/{dataset_meta.file}', 'r')
            self.stdout.write(f'GROUP COLUMNS: {list(adata.obs.keys())}')

        except BaseException as e:
            self.stdout.write(self.style.ERROR(str(e)))
