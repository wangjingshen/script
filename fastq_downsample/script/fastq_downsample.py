import pandas as pd
import subprocess
import sys
import argparse

from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

ROOT = Path(__file__).resolve().parent
#ROOT = os.path.dirname(os.path.abspath(__file__))
dev_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(dev_root))
from utils.utils import mkdir, logger, execute_cmd, timer, run_with_single_thread


@timer
def fq_subset(library_id, fq_path, spname, sub_G):
    mkdir(spname)

    if sub_G < 1000:
        num_reads = round(sub_G*10**9/150/2)
    else:
        num_reads = sub_num
    
    logger.info(f'[{spname}] subset {sub_G} G reads, {num_reads} reads.')
    cmd1 = f"seqtk sample -s 100 {fq_path}/{library_id}*R1*gz {num_reads} | gzip > {spname}/{library_id}_sub_{sub_G}G_R1.fq.gz"
    cmd2 = f"seqtk sample -s 100 {fq_path}/{library_id}*R2*gz {num_reads} | gzip > {spname}/{library_id}_sub_{sub_G}G_R2.fq.gz"
    try:
        execute_cmd(cmd1)
        execute_cmd(cmd2)
    except subprocess.CalledProcessError as e:
        print(f"Command failed with error: {e}")
        sys.exit(1)


def parse_mapfile(mapfile):
    df = pd.read_csv(mapfile, sep=r'\s+')
    columns = ['library_id', 'fq_path', 'spname', 'sub_G']    
    return [df[col] for col in columns]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mapfile', help='mapfile', required=True)
    # add version
    parser.add_argument('--version', action='version', version='1.0')

    args = parser.parse_args()
    library_id_list, fq_path_list, spname_list, sub_G_list = parse_mapfile(args.mapfile)

    with ProcessPoolExecutor(max_workers = 1) as executor:
        executor.map(fq_subset, library_id_list, fq_path_list, spname_list, sub_G_list)


if __name__ == '__main__':
    main()