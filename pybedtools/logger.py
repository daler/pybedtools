"""
Configures a logger to be used across all of pybedtools
"""
import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setLevel(logging.WARNING)
formatter = logging.Formatter("%(name)s [%(levelname)s]: %(message)s")
ch.setFormatter(formatter)
logger.addHandler(ch)
