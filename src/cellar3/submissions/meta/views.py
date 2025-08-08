import logging
import uuid
from typing import Dict

from rest_framework.response import Response
from rest_framework.views import APIView

from cellar3.submissions.proxy import ProxyDev
from src.cellar3.exceptions import UserException
from src.cellar3.models import Submission, SubmissionStatus, ConversionStatus
from src.cellar3.submissions.mailing import send_review_submission_email
from src.cellar3.tools import raise_if_not_dev, is_dev

logger = logging.getLogger('cellar.submit.meta')


class SubmitDatasetMetaView(APIView):
    META_MANDATORY_FIELDS = ["id", "file", "name", "tissue", "species", "cells", "public", "type", "info", "groups",
                             "categories"]

    def post(self, request, *args, **kwargs):
        """
        Submits the pre-submission inquiry:
        - Creates a submission record in the DB
        - Sends emails to managers/admins

        :param request:
        :param args:
        :param kwargs:
        :return:
        """
        if not is_dev():
            return ProxyDev.proxy_view(request)

        raise_if_not_dev()

        try:
            json_data = request.data
            meta = json_data.get('meta')
            notes = json_data.get('notes')

            if not meta:
                raise UserException('No meta JSON received.')

            self._qc(meta)

            submission_id = str(uuid.uuid4())
            self._create_submission(submission_id=submission_id, meta=meta)
            send_review_submission_email(submission_id=submission_id, notes=notes)

            return Response({"details": "Submission inquiry is being processed..."})
        except Exception as e:
            logger.exception(e)
            raise UserException(f'Could not process submission inquiry: {e}')

    @classmethod
    def _create_submission(cls, submission_id: str, meta: Dict):
        """
        :param submission_id:
        :param meta:
        :param notes:
        :return:
        """
        try:
            Submission(
                id=submission_id,
                maintainer_email=meta['maintainer']['email'],
                maintainer_name=meta['maintainer']['name'],
                status=SubmissionStatus.RECEIVED,
                conversion_status=ConversionStatus.INITIAL,
                data=meta
            ).save()
        except Exception as e:
            raise UserException(f'Could not process dataset: {e}')

    @classmethod
    def _qc(cls, meta: Dict):
        """
        @todo Make sure dataset ID doesnt already exist.
        :return:
        """
        for mandatory_field in cls.META_MANDATORY_FIELDS:
            if mandatory_field not in meta:
                raise UserException(f'Field {mandatory_field} is mandatory.')
