import os

from django.core.management.base import BaseCommand
from django.db.models import FileField

from src.cellar3.models import Submission, ConversionStatus


class Command(BaseCommand):
    help = 'Deletes the uploaded RDS file for a submission along with related records.'

    def add_arguments(self, parser):
        # Optional: Define command-line arguments here
        parser.add_argument('submission_id', type=str, help='ID of the submission')

    def handle(self, *args, **kwargs):
        try:
            submission_id = kwargs['submission_id']
            submission = Submission.objects.get(pk=submission_id)
            self._safe_delete_file(submission.rds)
            self._safe_delete_file(submission.h5ad)
            submission.error = None
            submission.rds = None
            submission.h5ad = None
            submission.conversion_status = ConversionStatus.INITIAL
            submission.save()
            self.stdout.write(self.style.SUCCESS(f'RDS data deleted for submission {submission_id}'))
        except Exception as e:
            self.stdout.write(self.style.ERROR(e))

    @classmethod
    def _safe_delete_file(cls, file_field: FileField):
        if file_field and os.path.exists(file_field.path):
            os.remove(file_field.path)
