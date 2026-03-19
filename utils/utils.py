import sys
import subprocess
import logging
import glob
import os
import time
import psutil
import time
#import functools
from functools import wraps
from contextlib import contextmanager
from pathlib import Path



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


def run_with_single_thread(command, **kwargs):
    '''
    Run the command in single-thread mode without affecting the main process.
    Enforce single-threading to avoid conflicts with OpenBLAS.
    '''
    logger.info(f"Executing: {command} with single thread")
    env = os.environ.copy()
    env.update({
        'OPENBLAS_NUM_THREADS': '1',
        'MKL_NUM_THREADS': '1',
        'OMP_NUM_THREADS': '1',
        'NUMEXPR_NUM_THREADS': '1',
    })

    try:
        subprocess.check_call(command, env=env, shell=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed: {command}, error: {e}")
        raise

    #return subprocess.run(cmd, env=env, **kwargs)


def timer(func):
    '''
    decorator: measure execution time
    '''
    @wraps(func)
    def wrapper(*args, **kwargs):
        logger.info(f"[{func.__name__}] start...")
        start = time.time()
        result = func(*args, **kwargs)
        elapsed = time.time() - start
        logger.info(f"[{func.__name__}] done. time: {elapsed:.2f}s ({elapsed/60:.2f}min) ({elapsed/60/60:.2f}h) ({elapsed/60/60/24:.2f}d)")
        return result
    return wrapper


# Context manager: temporarily switch the working directory
@contextmanager
def tmp_chdir(path: Path):
    origin = Path.cwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(origin)