from django.core.management import BaseCommand

from src.cellar3.exceptions import UserException
from src.cellar3.models import Submission
from src.cellar3.submissions.dataset.importer import ImporterService


class Command(BaseCommand):
    help = 'Imports meta from submission into Datasets table and copies the h5ad file to datasets folder.'

    def add_arguments(self, parser):
        # Optional: Define command-line arguments here
        parser.add_argument('submission_id', type=str, help='ID of the submission io import.')

    def handle(self, *args, **kwargs):
        try:
            submission_id = kwargs['submission_id']
            submission = Submission.objects.get(pk=submission_id)
            if not submission.h5ad.path:
                raise UserException('Could not import submission because it has no h5ad file.')
            ImporterService.import_submission_into_dev(submission_id=submission_id)
            self.stdout.write(self.style.SUCCESS(f'Successfully imported.'))
        except Exception as e:
            self.stdout.write(self.style.ERROR(str(e)))
