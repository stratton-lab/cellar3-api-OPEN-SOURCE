import logging
import time
from typing import Tuple

from django.conf import settings

from src.cellar3.exceptions import UserException

logger = logging.getLogger('cellar.timeit')


def _get_class_name(args: Tuple) -> str:
    if args:
        if hasattr(args[0], '__class__'):
            # This covers instances (methods) and types (class methods)
            class_name = args[0].__name__ if isinstance(args[0], type) else args[0].__class__.__name__
            class_name += '.'
        else:
            class_name = ''
    else:
        class_name = ''
    return class_name


def timeit(func):
    """
    A decorator that logs the execution time of the decorated function.
    """

    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)  # Call the actual function
        end_time = time.time()
        execution_time = end_time - start_time
        class_name = _get_class_name(args)

        logger.info(f"Function '{class_name}{func.__name__}' executed in {execution_time:.2f}s")
        return result

    return wrapper


def is_dev() -> bool:
    """
    Returns true if the server is deployed in the DEV server.
    :return:
    """
    return settings.DJANGO_ENV == "development"


def raise_if_not_dev():
    """
    Raises a UserException if function accessed outside a DEV server.
    :return:
    """
    if not is_dev():
        raise UserException('Forbidden. This feature is not available externally.')
