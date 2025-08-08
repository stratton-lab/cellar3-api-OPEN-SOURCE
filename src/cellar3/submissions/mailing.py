import logging

from django.conf import settings

from src.cellar3.email import send_emails
from src.cellar3.exceptions import UserException
from src.cellar3.models import Submission

logger = logging.getLogger('cellar.submit.mailing')


def send_review_submission_email(submission_id: str, notes: str):
    try:
        submission = Submission.objects.get(pk=submission_id)
        maintainer = submission.data.get('maintainer', {})

        submission_url = f'{settings.BACKEND_DEV_BASE_URL}/manage/cellar3/submission/{submission_id}'
        send_emails(
            emails=settings.SUBMIT_NOTIFY_EMAILS,
            subject='SingloCell - New dataset submitted',
            template='submission_review_email.html',
            context={
                'submission_id': submission_id,
                'maintainer_name': maintainer['name'],
                'maintainer_affiliation': maintainer['affiliation'],
                'maintainer_email': maintainer['email'],
                'submission_note': notes,
                'submission_url': submission_url,
            }
        )
    except Exception as e:
        raise UserException(f'Could not send review email: {e}')


def send_conversion_finished_email(submission_id: str):
    """
    Send email to maintainers that conversion is successful, with link to dataset on DEV.

    :param submission_id:
    :return:
    """
    try:
        submission = Submission.objects.get(pk=submission_id)
        dataset_id = submission.data['id']
        dataset_name = submission.data['name']
        link_to_dataset = f'{settings.FRONTEND_DEV_BASE_URL}display?dataset={dataset_id}'
        link_to_prod_dataset = f'{settings.FRONTEND_PROD_BASE_URL}display?dataset={dataset_id}'
        logger.info(f'Sending email with link to converted dataset: {link_to_dataset}')
        send_emails(
            emails=settings.SUBMIT_NOTIFY_EMAILS,
            subject='SingloCell - Submitted dataset Converted & Imported into DEV',
            template='submission_conversion_success_email.html',
            context={
                'submission_id': submission_id,
                'dataset_id': dataset_id,
                'dataset_name': dataset_name,
                'link_to_dataset': link_to_dataset,
                'link_to_prod_dataset': link_to_prod_dataset
            })
    except Exception as e:
        raise UserException(f'Could not send conversion finished email: {e}')


def send_conversion_failed_email(submission_id: str):
    """

    :param submission_id:
    :return:
    """
    try:
        submission = Submission.objects.get(pk=submission_id)
        dataset_id = submission.data['id']
        dataset_name = submission.data['name']
        logger.info(f'Sending email that conversion failed: {submission.error}')
        send_emails(
            emails=settings.SUBMIT_NOTIFY_EMAILS,
            subject='SingloCell - Submitted dataset could NOT be converted.',
            template='submission_conversion_failed_email.html',
            context={
                'submission_id': submission_id,
                'dataset_id': dataset_id,
                'dataset_name': dataset_name,
                'conversion_error': submission.error
            }
        )
    except Exception as e:
        raise UserException(f'Could not send conversion failed email: {e}')
