import logging
import sys


def setup_logger():
    logger = logging.getLogger("diamond_prep")
    logger.setLevel(logging.INFO)
    handler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger


__all__ = ["setup_logger"]
