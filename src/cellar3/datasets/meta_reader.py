import json
from typing import Dict, List

import yaml
from django.conf import settings

from src.cellar3.exceptions import UserException
from src.cellar3.models import DatasetMeta
from src.cellar3.serializers import DatasetInfoSerializer


class AbstractMetaReader:
    @classmethod
    def get_dataset_info(cls, dataset_id: str) -> Dict:
        pass

    @classmethod
    def get_public_datasets(cls) -> Dict[str, any]:
        pass

    @classmethod
    def get_multiple_dataset_info(cls, dataset_ids: List[str]) -> List[dict]:
        pass


class JSONMetaReader(AbstractMetaReader):
    """
    Loads datasets meta data from a JSON file.
    """

    @classmethod
    def _load_datasets_info(cls) -> List[Dict]:
        """
        Returns all datasets, both public and private.

        :return:
        """
        datasets_info = json.load(open(settings.DATASETS_META_FILE, 'r'))

        for dataset_info in datasets_info:
            dataset_info.pop('maintainer', None)

        return datasets_info

    @classmethod
    def get_dataset_info(cls, dataset_id: str) -> Dict:
        """
        Returns metadata for a dataset
        :param dataset_id:
        :return:
        """
        dataset = next((dataset for dataset in cls._load_datasets_info() if dataset['id'] == dataset_id), None)
        if not dataset:
            raise UserException(f'No datasets available for id {dataset_id}.')
        return dataset

    @classmethod
    def get_public_datasets(cls) -> Dict[str, any]:
        all_datasets = cls._load_datasets_info()
        available_datasets = [dataset for dataset in all_datasets if cls._should_display_dataset(dataset)]
        return {
            'datasets': available_datasets,
            'categories': yaml.safe_load(open('categories.yaml', 'r'))
        }

    @classmethod
    def _should_display_dataset(cls, meta: Dict):
        """
        We display the dataset only if it's public or setting sallow us to display private datasets.
        :param meta:
        :return:
        """
        return meta.get('public', False) or settings.SHOW_PRIVATE_DATASETS


class SQLiteMetaReader(AbstractMetaReader):
    """
    Loads datasets meta data from a SQLite database.
    """

    @classmethod
    def get_dataset_info(cls, dataset_id: str) -> Dict:
        try:
            dataset_meta = DatasetMeta.objects.get(pk=dataset_id)
            serializer = DatasetInfoSerializer(dataset_meta)
            return serializer.data
        except DatasetMeta.DoesNotExist:
            raise UserException(f'No datasets available for id {dataset_id}.')
        except BaseException:
            raise UserException(f'Could not load dataset for id {dataset_id}.')

    @classmethod
    def get_public_datasets(cls) -> List[Dict]:
        objects = DatasetMeta.objects
        queryset = objects if settings.SHOW_PRIVATE_DATASETS else objects.filter(public=True)
        serializer = DatasetInfoSerializer(queryset, many=True)
        return serializer.data

    @classmethod
    def get_multiple_dataset_info(cls, dataset_ids: List[str]) -> List[Dict]:
        dataset_meta_qs = DatasetMeta.objects.filter(pk__in=dataset_ids)
        serializer = DatasetInfoSerializer(dataset_meta_qs, many=True)
        return serializer.data


# Allows to easily switch between JSON file reader and Sqlite reader
META_READER = SQLiteMetaReader
