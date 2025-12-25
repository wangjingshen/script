import sys
import subprocess
import logging
import glob
import os


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)


def execute_cmd(command) -> None:
    logger.info(f"Executing: {command}")
    try:
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed: {command}, error: {e}")
        raise


def find_file(pattern: str) -> str:
    """
    Args: pattern of file 
    Returns: file path
    """
    files = glob.glob(pattern)
    if not files:
        raise PathError(f"file not find: {pattern}")
    return files[0]


def mkdir(dir) -> None:
    try:
        os.makedirs(dir, exist_ok=True)
    except Exception as e:
        logger.error(f"Failed to create directory {dir}: {e}")