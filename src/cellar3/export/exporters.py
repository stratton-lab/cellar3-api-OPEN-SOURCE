import abc
import logging
from io import StringIO, BytesIO
from typing import Dict

from pandas import DataFrame

logger = logging.getLogger('cellar.export.exporters')


class Exporter:
    extension = None  # ex: csv
    content_type = None  # ex: application/vnd.ms-excel

    @abc.abstractmethod
    def get_content(self, params: Dict[str, any]) -> str | bytes:
        pass


class CSVExporter(Exporter):
    extension = 'csv'
    content_type = 'text/csv'

    def get_content(self, df: DataFrame) -> str:
        buffer = StringIO()
        df.to_csv(buffer, index=False)
        return buffer.getvalue()


class ExcelExporter(Exporter):
    extension = 'xlsx'
    content_type = 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'

    def get_content(self, df: DataFrame) -> str | bytes:
        buffer = BytesIO()
        df.to_excel(buffer, index=False)
        return buffer.getvalue()


class JSONExporter(Exporter):
    extension = 'json'
    content_type = 'application/json'

    def get_content(self, df: DataFrame) -> str | bytes:
        json_data = df.to_json(orient='records')
        return json_data
