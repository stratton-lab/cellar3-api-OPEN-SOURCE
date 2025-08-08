from rest_framework import status
from rest_framework.exceptions import APIException


class UserException(APIException):
    status_code = status.HTTP_400_BAD_REQUEST
    default_detail = 'Could not process user request.'
    default_code = 'invalid'


class ExternalResourceException(APIException):
    status_code = status.HTTP_400_BAD_REQUEST
    default_detail = 'Could not retrieve data.'
    default_code = 'invalid'
