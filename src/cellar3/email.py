import logging
import smtplib
from email.message import EmailMessage
from typing import List, Dict

from django.conf import settings
from django.template.loader import render_to_string
from premailer import transform as html_transform

from src.cellar3.exceptions import UserException

logger = logging.getLogger('cellar.emails')


def send_emails(emails: List[str], subject: str, template: str, context: Dict):
    if not settings.DEBUG:
        for email in emails:
            send_email(email=email, subject=subject, template=template, context=context)
    else:
        logger.info('Skipping sending emails because running in local')


def send_email(email: str, subject: str, template: str, context: Dict):
    """
    Email server: EMAIL_SERVER_URL

    :param context: Dict with parameters to fill the template.
    :param template: Name of template file
    :param email:
    :param subject:
    :return:
    """
    try:
        msg = EmailMessage()
        msg['Subject'] = subject
        msg['From'] = 'EMAIL_SERVER_URL'
        msg['To'] = email

        html_content = render_to_string(template, context)
        html_content_inlined = html_transform(html_content)
        msg.set_content(html_content_inlined, subtype='html')

        server = smtplib.SMTP(settings.EMAIL_SERVER)
        server.send_message(msg)
        server.quit()
    except Exception as e:
        logger.error(f'Could not send email: {e}')
        raise UserException(f'Could not send email.')
