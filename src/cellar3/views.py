import logging

from django.conf import settings
from rest_framework.response import Response
from rest_framework.views import APIView

logger = logging.getLogger('cellar.view')


class IndexView(APIView):
    def get(self, request):
        return Response({
            'status': 'OK',
            'environment': settings.DJANGO_ENV
        })
