import logging
import os

import requests
from django.conf import settings
from django.core.files.base import ContentFile

from src.cellar3.exceptions import UserException
from src.cellar3.models import Submission, ConversionStatus

logger = logging.getLogger('cellar.submit.convert')


class ConversionService:

    @classmethod
    def convert(cls, submission_id: str):
        submission = Submission.objects.get(pk=submission_id)
        try:
            rds_path = submission.rds.path
            cls._set_status(submission, ConversionStatus.RUNNING)
            files = {'file': open(rds_path, 'rb')}
            logger.info('Calling conversion service...')
            r = requests.post(f'{settings.CONVERSION_SERVICE_URL}/convert/', files=files)
            if not r.status_code == requests.codes.ok:
                raise UserException(f'Conversion failed: {r.json()["detail"]}')
            file_content = ContentFile(r.content, cls._get_h5ad_name(rds_path))
            logger.info('Conversion successful.')
            cls._set_status_success(submission, file_content)
        except requests.exceptions.ConnectionError as e:
            logger.exception(e)
            cls._set_status_error(submission, e)
            raise UserException('Conversion service is not available.')
        except Exception as e:
            logger.exception(e)
            cls._set_status_error(submission, e)
            raise UserException(f'Conversion Service returned an error: {e}')

    @classmethod
    def _get_h5ad_name(cls, rds_path: str) -> str:
        base_path, original_ext = os.path.splitext(os.path.basename(rds_path))
        return f"{base_path}.h5ad"

    @classmethod
    def _set_status(cls, submission: Submission, status: ConversionStatus):
        submission.conversion_status = status
        submission.save()

    @classmethod
    def _set_status_error(cls, submission: Submission, e: Exception):
        submission.conversion_status = ConversionStatus.ERROR
        submission.error = str(e)
        submission.save()

    @classmethod
    def _set_status_success(cls, submission: Submission, file_content: ContentFile):
        submission.conversion_status = ConversionStatus.SUCCESS
        submission.h5ad = file_content
        submission.error = None
        submission.save()
        logger.info('H5AD file saved.')
