import abc
import logging
from typing import Dict

from django.http import HttpResponse
from pandas import DataFrame
from rest_framework.views import APIView

from src.cellar3.exceptions import UserException
from src.cellar3.export.exporters import CSVExporter, ExcelExporter, JSONExporter

logger = logging.getLogger('cellar.export.views')


class ExportView(APIView):
    exporters = {
        'csv': CSVExporter(),
        'excel': ExcelExporter(),
        'json': JSONExporter()
    }

    def post(self, request, dataset_id: str, export_format: str):
        return self._export(dataset_id, export_format, request.data)

    def _export(self, dataset_id: str, export_format: str, params: Dict):
        exporter = self.exporters.get(export_format)
        if not exporter:
            raise UserException(f'Unsupported export format. Supported formats: {list(self.exporters.keys())}')
        response = HttpResponse(content_type=exporter.content_type)
        response['Content-Disposition'] = f'attachment; filename="{self.get_file_name(params)}.{exporter.extension}"'
        df = self.get_df(dataset_id, params)
        response.write(exporter.get_content(df))
        return response

    @abc.abstractmethod
    def get_file_name(self, params: Dict[str, any]) -> str:
        pass

    @abc.abstractmethod
    def get_df(self, dataset_id: str, params: Dict[str, any]) -> DataFrame:
        pass
