import logging
import os
from concurrent.futures import ThreadPoolExecutor
from typing import List, Dict

from django.conf import settings
from django.core.files.uploadedfile import UploadedFile
from rest_framework.request import Request
from rest_framework.response import Response
from rest_framework.views import APIView

from cellar3.submissions.proxy import ProxyDev
from src.cellar3.exceptions import UserException
from src.cellar3.models import Submission, DatasetMeta, ConversionStatus
from src.cellar3.submissions.dataset.conversion import ConversionService
from src.cellar3.submissions.dataset.importer import ImporterService
from src.cellar3.submissions.mailing import send_conversion_finished_email, send_conversion_failed_email
from src.cellar3.tools import raise_if_not_dev, is_dev

logger = logging.getLogger('cellar.submit.dataset')


class SubmitDatasetFileView(APIView):
    _executor = ThreadPoolExecutor()  # Persistent executor shared across requests

    @classmethod
    def _get_submission(cls, submission_id: str) -> Submission:
        """
        Checks that the user can retrieve the submission, and returns the submission record from DB.

        :param submission_id:
        :return:
        """
        try:
            logger.info(f'Retrieving submission {submission_id}')
            submission = Submission.objects.get(pk=submission_id)
            if submission.rds:
                raise UserException('An RDS file has already been uploaded for this submission')

            return submission

        except Submission.DoesNotExist:
            raise UserException('There are no existing submissions with this id. '
                                'Please make sure to fill-out the pre-submission inquiry form first.')
        except Exception as e:
            raise UserException(f'Can not access submission data: {e}')

    def get(self, request, submission_id: str):
        """
        Allows to retrieve basic info about a submission.

        :param request:
        :param submission_id:
        :return:
        """
        if not is_dev():
            return ProxyDev.proxy_view(request)  # Returns result from DEV server.

        raise_if_not_dev()
        submission = self._get_submission(submission_id=submission_id)
        self._check_unique_dataset(submission)
        try:
            return Response({
                'id': submission.id,
                'dataset_id': submission.data['id'],
                'dataset_name': submission.data['name'],
                'spatial_embeddings': self.get_spatial_embeddings(submission)
            })
        except Exception:
            raise UserException('Could not retrieve submission data.')

    def post(self, request, submission_id: str):
        """
        Allows to upload an RDS dataset file.
        @todo Accept additional file types: rds, h5ad, cloupe

        :param request:
        :param submission_id:
        :return:
        """
        if not is_dev():
            return ProxyDev.proxy_view(request)

        raise_if_not_dev()
        submission = self._get_submission(submission_id=submission_id)
        rds_file = request.FILES.get('file')

        if rds_file:
            self._check_uploaded_file(rds_file)
            self._check_unique_dataset(submission)
            try:
                # Saving RDS file
                dataset_id = submission.data['id']
                rds_file.name = f'{dataset_id}.rds'
                submission.rds = rds_file
                submission.error = None
                submission.conversion_status = ConversionStatus.INITIAL
                self.save_preview_image(submission, request)
                self.save_spatial_images(submission, request)
                submission.save()
                logger.info('RDS file saved.')

                # Converting to H5AD (Running in Thread)
                self.convert_and_import_in_background(submission_id=submission_id)

                return Response({"details": "Your RDS file has been successfully uploaded. "
                                            "We will contact you once your dataset is available on the website."})
            except Exception as e:
                submission.error = str(e)
                submission.save()
                send_conversion_failed_email(submission_id)
                raise UserException(f'Could not process dataset: {e}')
        else:
            raise UserException("No file uploaded.")

    def convert_and_import_in_background(self, submission_id: str):
        self._executor.submit(SubmitDatasetFileView.convert_and_import, submission_id=submission_id)

    @staticmethod
    def convert_and_import(submission_id: str):
        try:
            ConversionService.convert(submission_id=submission_id)
            ImporterService.import_submission_into_dev(submission_id=submission_id)
            send_conversion_finished_email(submission_id=submission_id)
        except Exception as e:
            logger.exception(e)
            submission = Submission.objects.get(pk=submission_id)
            submission.error = str(e)
            submission.save()
            send_conversion_failed_email(submission_id)

    @staticmethod
    def _check_uploaded_file(rds_file: UploadedFile):
        """
        @todo Check for a list of supported file types
        :param rds_file:
        :return:
        """
        base_path, ext = os.path.splitext(rds_file.name)
        if not ext == '.rds':
            raise UserException('Only RDS files are supported.')

    @staticmethod
    def _check_unique_dataset(submission: Submission):
        """
        @todo Check that dataset_id doesn already exist in DatasetMeta table.
        :param submission:
        :return:
        """
        dataset_id = submission.data['id']
        if DatasetMeta.objects.filter(id=dataset_id).exists():
            raise UserException('A dataset with this id already exists.')

    @classmethod
    def save_preview_image(cls, submission: Submission, request: Request):
        try:
            preview_image_file = request.FILES.get('previewImage')
            if preview_image_file:
                ext = os.path.splitext(preview_image_file.name)[1].lower()  # e.g. ".png"
                image_file_name = f"{submission.data['id']}{ext}"
                file_path = os.path.join(settings.MEDIA_ROOT, "datasets", image_file_name)
                cls.save_uploaded_file_to_disk(preview_image_file, file_path)
                submission.data['image'] = image_file_name
                logger.info('Preview image saved.')
            else:
                logger.warning('No preview image submitted by user.')
        except Exception as e:
            logger.warning('Could not save preview image.')
            logger.exception(e)

    @classmethod
    def save_spatial_images(cls, submission: Submission, request: Request):
        try:
            spatial_images = request.FILES.getlist('spatialImages')
            if spatial_images:
                spatial_embeddings = cls.get_spatial_embeddings(submission)
                cls.qc_spatial_embeddings(spatial_images, spatial_embeddings)
                for i, spatial_image in enumerate(spatial_images):
                    spatial_embeddings[i]['image'] = spatial_image.name

                    # Saving file to disk
                    folder_path = os.path.join(settings.MEDIA_ROOT, 'spatial', submission.data['id'])
                    os.makedirs(folder_path, exist_ok=True)
                    file_path = os.path.join(str(folder_path), spatial_image.name)
                    cls.save_uploaded_file_to_disk(spatial_image, file_path)
                    logger.info(f'Spatial image {spatial_image.name} saved.')
            else:
                logger.warning('No spatial images submitted by user.')
        except Exception as e:
            logger.warning('Could not save spatial images.')
            logger.exception(e)

    @classmethod
    def get_spatial_embeddings(cls, submission: Submission) -> List[Dict]:
        return [e for e in submission.data.get('embeddings', []) if e.get('key', '').startswith('spatial_')]

    @classmethod
    def qc_spatial_embeddings(cls, spatial_images: List[UploadedFile], embeddings: List[Dict]):
        if len(spatial_images) != len(embeddings):
            raise UserException(
                f'Number of spatial images ({len(spatial_images)}) does not match number of embeddings ({len(embeddings)}).')

    @classmethod
    def save_uploaded_file_to_disk(cls, file: UploadedFile, file_path: str):

        with open(file_path, 'wb+') as f:
            for chunk in file.chunks():
                f.write(chunk)
