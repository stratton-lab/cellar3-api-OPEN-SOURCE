from django.core.management import BaseCommand

from src.cellar3.models import DatasetMeta

ORANGE = "\033[33m"
GREEN = "\033[32m"
BLUE = "\033[36m"
RESET = "\033[0m"


class Command(BaseCommand):
    help = 'Lists all datasets with their IDs and names.'

    def handle(self, *args, **kwargs):
        try:
            datasets = DatasetMeta.objects.all()

            if not datasets.exists():
                self.stdout.write(self.style.WARNING('No datasets found.'))
                return

            self.stdout.write(self.style.SUCCESS(f'\nListing {datasets.count()} datasets:\n'))
            max_id_length = max(len(str(dataset.id)) for dataset in datasets)
            id_width = max(max_id_length, 10)  # Set a minimum width of 10 characters

            for dataset in datasets:
                dataset_id = str(dataset.id).ljust(id_width)
                dataset_name = dataset.name
                status = f'{ORANGE}PUBLIC{RESET}' if dataset.public else f'{GREEN}private{RESET}'

                self.stdout.write(f'{BLUE}{dataset_id}{RESET}\t{status}\t{dataset_name}')

        except Exception as e:
            self.stdout.write(self.style.ERROR(f'Error: {str(e)}'))
