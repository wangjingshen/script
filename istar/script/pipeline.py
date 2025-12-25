import sys
import pandas as pd
import argparse
import os
import glob
import subprocess
import logging
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

ROOT = Path(__file__).resolve().parent
#ROOT = os.path.dirname(os.path.abspath(__file__))

CONFIG = {
    "DEVICE": "cpu",  # or "cuda"
    "RADIUS_RAW": 20,
    "PixelSizeRaw": 0.25,
    "PixelSize": 0.5,
    "N_HVG": 1000,
    "EPOCHS": 400,
    "FilterSize": 8,
    "MinClusterSize": 20,
    "N_CLUSTERS": 10
}


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


class IStar:
    def __init__(self, dir: Path, image: Path, pos: Path, spname: str, step: str):
        filtered_dir = os.path.join(dir, "outs", "filtered")
        if not os.path.exists(filtered_dir):
            raise FileNotFoundError(f"[IStar] Required directory not found: {filtered_dir}")
        self.mtx = filtered_dir

        if not os.path.isfile(image):
            raise FileNotFoundError(f"[IStar] Image file not found: {image}")
        if not os.path.isfile(pos):
            raise FileNotFoundError(f"[IStar] Position file not found: {pos}")
        self.image = image
        self.pos = pos
        self.spname = spname.rstrip("/") + "/"  #  spname+"/"  #istar prefix
        self.step = step.strip().split(',')
        
        try:
            os.makedirs(self.spname, exist_ok=True)
        except Exception as e:
            logger.error(f"Failed to create directory {self.spname}: {e}")

    def input(self) -> None:
        '''
        generate_input
        '''
        logger.info(f"[{self.spname}] Running input step.")
        execute_cmd((f'/SGRNJ/Public/Software/conda_env/r4.1_env/bin/Rscript {ROOT}/gen_cnts_locs.R '
                f'--mtx {self.mtx} '
                f'--pos {self.pos} '
                f'--spname {self.spname} '
                ))
        if not os.path.exists("checkpoints"):
            execute_cmd(f'cp -r {ROOT}/../data/checkpoints/ .')

        with open(f"{self.spname}/radius-raw.txt", "w") as f:
            f.write(str(CONFIG["RADIUS_RAW"]))
        with open(f"{self.spname}/pixel-size-raw.txt", "w") as f:
            f.write(str(CONFIG["PixelSizeRaw"]))
        with open(f"{self.spname}/pixel-size.txt", "w") as f:
            f.write(str(CONFIG["PixelSize"]))

        execute_cmd(f'cp {self.image} {self.spname}/he-raw.png')
        logger.info(f"[{self.spname}] input done.")

    def preprocess(self) -> None:
        '''
        preprocess_image
        '''
        logger.info(f"[{self.spname}] Running preprocess step.")
        execute_cmd(f'python {ROOT}/istar/rescale.py {self.spname}/ --image')
        execute_cmd(f'python {ROOT}/istar/preprocess.py {self.spname} --image')
        # extract histology features
        execute_cmd(f'python {ROOT}/istar/extract_features.py {self.spname} --device={CONFIG["DEVICE"]}')
        # auto detect tissue mask
        execute_cmd(f'python {ROOT}/istar/get_mask.py {self.spname}/embeddings-hist.pickle {self.spname}/mask-small.png')
        # select most highly variable genes to predict
        execute_cmd(f'python {ROOT}/istar/select_genes.py --n-top={CONFIG["N_HVG"]} {self.spname}/cnts.tsv {self.spname}/gene-names.txt')
        # rescale coordinates and spot radius
        execute_cmd(f'python {ROOT}/istar/rescale.py {self.spname} --locs --radius')
        logger.info(f"[{self.spname}] preprocess done.")

    def impute(self) -> None:
        '''
        predict super-resolution gene expression
        '''
        logger.info(f"[{self.spname}] Running impute step.")
        # train gene expression prediction model and predict at super-resolution
        execute_cmd(f'python {ROOT}/istar/impute.py {self.spname} --epochs={CONFIG["EPOCHS"]} --device={CONFIG["DEVICE"]}')
        # visualize imputed gene expression
        execute_cmd(f'python {ROOT}/istar/plot_imputed.py {self.spname}')
        logger.info(f"[{self.spname}] impute done.")

    def cluster(self) -> None:
        logger.info(f"[{self.spname}] Running cluster step.")
        # segment image by gene features
        execute_cmd(f'python {ROOT}/istar/cluster.py --filter-size={CONFIG["FilterSize"]} --min-cluster-size={CONFIG["MinClusterSize"]} --n-clusters={CONFIG["N_CLUSTERS"]} --mask={self.spname}/mask-small.png {self.spname}/embeddings-gene.pickle {self.spname}/clusters-gene/')
        # differential analysis by clusters
        execute_cmd(f'python {ROOT}/istar/aggregate_imputed.py {self.spname}')
        execute_cmd(f'python {ROOT}/istar/reorganize_imputed.py {self.spname}')
        execute_cmd(f'python {ROOT}/istar/differential.py {self.spname}')
        # visualize spot-level gene expression data
        execute_cmd(f'python {ROOT}/istar/plot_spots.py {self.spname}')
        logger.info(f"[{self.spname}] cluster done.")

    def run(self) -> None:
        step_order = ['input','preprocess','impute','cluster']
        for step in step_order:
            if step in self.step:
                getattr(self, step)()


def parse_mapfile(mapfile, step):
    df_mapfile = pd.read_csv(mapfile, sep=r'\s+')
    df_mapfile['step'] = step
    dir_list = df_mapfile['dir']
    image_list = df_mapfile['image']
    pos_list = df_mapfile['pos']
    spname_list = df_mapfile['spname']
    step_list = df_mapfile['step']
    return dir_list, image_list, pos_list, spname_list, step_list

def run_single(dir, image, pos, spname, step):
    try:
        runner = IStar(dir, image, pos, spname, step)
        runner.run()
        logger.info(f'Completed: {spname}')
    except Exception as e:
        logger.error(f'JOB FAILED: {spname}',
              file=sys.stderr)
        traceback.print_exc()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mapfile', help='tsv:dir<>image<>pos<>spname', required=True)
    parser.add_argument('--step', default='input,preprocess,impute,cluster', help='comma-separated step')
    parser.add_argument('--threads', type=int, default=1, help='thread pool size')
    args = parser.parse_args()

    dir_list, image_list, pos_list, spname_list, step_list = parse_mapfile(args.mapfile, args.step)
    logger.info(f"Starting pipeline with {len(dir_list)} samples and {args.threads} thread(s).")
    with ThreadPoolExecutor(max_workers = args.threads) as executor:
        executor.map(run_single, dir_list, image_list, pos_list, spname_list, step_list)


if __name__ == '__main__':
    main()