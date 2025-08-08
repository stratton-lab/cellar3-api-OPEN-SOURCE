import json
from functools import lru_cache
from typing import List, Tuple, Dict

from pandas import DataFrame

from src.cellar3.analysis.differential.dge import DifferentialExpression
from src.cellar3.analysis.differential.subset import ExpressionMatrix
from src.cellar3.datasets.dataset import Dataset
from src.cellar3.receipt import Receipt


def _unstringify_ints(list_str: str):
    return [int(item) for item in (list_str or '').split(',')]


def _stringify_ints(items: List[int]):
    return ','.join([str(item) for item in items])


@lru_cache(maxsize=20)
def _cached_dge(dataset_id: str, target_str: str, background_str: str, diff_p_val_threshold: float,
                fold_change_threshold: float, meta_str: str) -> Tuple[Dataset, Dict, DataFrame]:
    target = _unstringify_ints(target_str)
    background = _unstringify_ints(background_str)
    meta = json.loads(meta_str)

    dataset = Dataset(dataset_id)
    Receipt.add_dataset_load_step(meta, dataset.meta)

    adata = ExpressionMatrix.extract(dataset, target=target, background=background, meta=meta)
    return dataset, meta, DifferentialExpression(pval_thr=diff_p_val_threshold,
                                                 fold_change_threshold=fold_change_threshold) \
        .process(adata, target=target, background=background, meta=meta)


def get_cached_dge(dataset_id: str, target: List[int], background: List[int], diff_p_val_threshold: float,
                   fold_change_threshold: float, meta: Dict) -> Tuple[Dataset, Dict, DataFrame]:
    target_str = _stringify_ints(target)
    background_str = _stringify_ints(background)
    meta_str = json.dumps(meta)
    return _cached_dge(dataset_id=dataset_id, target_str=target_str, background_str=background_str,
                       diff_p_val_threshold=diff_p_val_threshold, fold_change_threshold=fold_change_threshold,
                       meta_str=meta_str)
