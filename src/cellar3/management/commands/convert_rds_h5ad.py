from django.core.management import BaseCommand

from src.cellar3.submissions.dataset.conversion import ConversionService


class Command(BaseCommand):
    help = 'Converts the RDS file in this submission into an H%AD file. Both files are saved on disk.'

    def add_arguments(self, parser):
        # Optional: Define command-line arguments here
        parser.add_argument('submission_id', type=str, help='ID of the submission with the RDS file to convert.')

    def handle(self, *args, **kwargs):
        try:
            submission_id = kwargs['submission_id']
            ConversionService.convert(submission_id=submission_id)

            self.stdout.write(self.style.SUCCESS(f'Successfully converted.'))
        except Exception as e:
            self.stdout.write(self.style.ERROR(e))
