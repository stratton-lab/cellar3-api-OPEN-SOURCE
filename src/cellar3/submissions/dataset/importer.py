import logging
import shutil

from django.conf import settings

from src.cellar3.exceptions import UserException
from src.cellar3.models import Submission, DatasetMeta

logger = logging.getLogger('cellar.submit.import')


class ImporterService:
    """
    Imports a submission into the DEV database.
    """

    @staticmethod
    def _import_meta_into_dev(submission_id: str):
        """
        Import data json from submission object into datasets table.

        :return:
        """
        try:
            logger.info('Importing meta data from submission into the DEV database...')
            submission = Submission.objects.get(pk=submission_id)
            DatasetMeta.objects.create(**submission.data)
        except Exception as e:
            logger.exception(e)
            raise UserException(f'Could not import data into datasets: {e}')

    @staticmethod
    def _import_h5ad_into_dev(submission_id: str):
        """
        Copy h5ad file from submissions folder to datasets folder.
        @todo Make sure django has writing rights to datasets folder
        :param submission_id:
        :return:
        """
        try:
            logger.info('Copying h5ad file from submissions folder into the DEV datasets folder...')
            submission = Submission.objects.get(pk=submission_id)
            shutil.copy(submission.h5ad.path, settings.DATASETS_FOLDER)
        except Exception as e:
            logger.exception(e)
            raise UserException(f'Could not copy h5ad file: {e}')

    @classmethod
    def import_submission_into_dev(cls, submission_id: str):
        cls._import_meta_into_dev(submission_id)
        cls._import_h5ad_into_dev(submission_id)
