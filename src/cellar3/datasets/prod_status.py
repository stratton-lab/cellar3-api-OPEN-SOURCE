import logging
from typing import Dict, List

import requests
from django.conf import settings

from cellar3.exceptions import UserException
from cellar3.tools import is_dev

logger = logging.getLogger('cellar.dataset.status.prod')

STATUS_DEV_ONLY = 'DEV_ONLY'
STATUS_PROD_PUBLIC = 'PROD_PUBLIC'
STATUS_PROD_PRIVATE = 'PROD_PRIVATE'

def get_datasets_prod_status(datasets: List[Dict]) -> Dict[str, str]:
    try:
        if not is_dev():
            return {}
        dataset_ids = [dataset['id'] for dataset in datasets]
        res = requests.post(f'{settings.BACKEND_PROD_BASE_URL}dataset/status/public', json=dataset_ids)
        res.raise_for_status()
        prod_datasets_public_status = res.json()
        return {d_id: _get_prod_status(d_id, prod_datasets_public_status) for d_id in dataset_ids}
    except Exception as e:
        logger.exception(e)
        raise UserException(f'Could not retrieve prod status of datasets: {e}')


def _get_prod_status(dataset_id, prod_statuses: Dict[str, bool]):
    if dataset_id not in prod_statuses:
        return STATUS_DEV_ONLY
    return STATUS_PROD_PUBLIC if prod_statuses[dataset_id] else STATUS_PROD_PRIVATE


def assign_prod_status(datasets: List[Dict]):
    if is_dev():
        try:
            datasets_prod_status = get_datasets_prod_status(datasets)
            for dataset in datasets:
                dataset['prodStatus'] = datasets_prod_status.get(dataset['id'])
        except Exception as e:
            logger.exception(e)
            raise UserException(f'Could not assign prod status to datasets: {e}')
