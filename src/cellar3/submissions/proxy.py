import logging

import requests
from django.conf import settings
from django.http import HttpResponse

from cellar3.exceptions import UserException

logger = logging.getLogger('cellar.submit.proxy')


class ProxyDev:
    """
    Proxies specific requests from PROD server to DEV server.
    For example dataset submission requests.
    """

    @classmethod
    def proxy_view(cls, request):
        try:
            target_url = f"{settings.BACKEND_DEV_BASE_URL}{request.get_full_path().replace(settings.URL_PREFIX, '')}"
            logger.info(f'Proxying {request.method} request to {target_url}')

            # Support for case where single key has multiple values (files)
            files = [
                (field, (file.name, file, file.content_type))
                for field in request.FILES
                for file in request.FILES.getlist(field)
            ] if request.FILES else None

            response = requests.request(
                method=request.method,
                url=target_url,
                json=request.data if request.content_type == 'application/json' and not request.FILES else None,
                files=files,
                params=request.GET.dict() if request.method == 'GET' else None
            )

            return HttpResponse(
                content=response.content,
                status=response.status_code,
                content_type=response.headers.get('Content-Type', 'application/json')
            )

        except BaseException as e:
            logger.exception(e)
            raise UserException(f'Could not proxy request: {e}')
