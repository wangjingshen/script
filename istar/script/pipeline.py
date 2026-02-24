import os
os.environ['OMP_NUM_THREADS'] = '32'
os.environ['MKL_NUM_THREADS'] = '32'
#os.environ['OPENBLAS_NUM_THREADS'] = '1'
#os.environ['NUMEXPR_NUM_THREADS'] = '1'
import sys
import pandas as pd
import argparse
import glob
import subprocess
import logging
import json
import psutil

from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

ROOT = Path(__file__).resolve().parent
#ROOT = os.path.dirname(os.path.abspath(__file__))
dev_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(dev_root))
from utils.utils import mkdir, logger, execute_cmd, timer, run_with_single_thread


CONFIG = {
    "DEVICE": "cpu",  # or "cuda"
    "PixelSize": 0.5,
    "N_HVG": 1000,
    "EPOCHS": 400,
    "FilterSize": 8,
    "MinClusterSize": 20,
    "N_CLUSTERS": 10
}



class IStar:
    def __init__(self, dir: Path, image: Path, spname: str, swap_pos: str, step: str):
        self.dir = dir
        filtered_dir = os.path.join(dir, "outs", "filtered")
        pos = os.path.join(dir, "outs", "spatial", "positions_list.csv")

        if not os.path.exists(filtered_dir):
            raise FileNotFoundError(f"[IStar] Required directory not found: {filtered_dir}")
        if not os.path.isfile(image):
            raise FileNotFoundError(f"[IStar] Image file not found: {image}")
        if not os.path.isfile(pos):
            raise FileNotFoundError(f"[IStar] Position file not found: {pos}")
        
        self.mtx = filtered_dir
        self.image = image
        self.pos = pos

        self.swap_pos = swap_pos
        self.spname = spname.rstrip("/") + "/"  #  spname+"/"  #istar prefix
        self.step = step.strip().split(',')


    @timer
    def input(self) -> None:
        '''
        generate_input
        '''
        logger.info(f"[{self.spname}] Running input step.")
        mkdir(self.spname)
        execute_cmd((f'/SGRNJ/Public/Software/conda_env/r4.1_env/bin/Rscript {ROOT}/gen_cnts_locs.R '
                f'--mtx {self.mtx} '
                f'--pos {self.pos} '
                f'--spname {self.spname} '
                f'--swap_pos {self.swap_pos}'
                ))
        if not os.path.exists("checkpoints"):
            execute_cmd(f'cp -r {ROOT}/../data/checkpoints/ .')

        with open(f'{self.dir}/outs/spatial/scalefactors_json.json') as f:
            sf = json.load(f)
        PixelSizeRaw = 8000/2000 * float(sf['tissue_hires_scalef'])  #https://github.com/daviddaiweizhang/istar
        RADIUS_RAW = 0.5 * float(sf['spot_diameter_fullres'])

        with open(f"{self.spname}/radius-raw.txt", "w") as f:
            f.write(str(RADIUS_RAW))
        with open(f"{self.spname}/pixel-size-raw.txt", "w") as f:
            f.write(str(PixelSizeRaw))
        with open(f"{self.spname}/pixel-size.txt", "w") as f:
            f.write(str(CONFIG["PixelSize"]))

        execute_cmd(f'cp {self.image} {self.spname}/he-raw.png')
        logger.info(f"[{self.spname}] input done.")

    @timer
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

    @timer
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

    @timer
    def cluster(self) -> None:
        logger.info(f"[{self.spname}] Running cluster step.")
        # segment image by gene features
        #execute_cmd(f'python {ROOT}/istar/cluster.py --filter-size={CONFIG["FilterSize"]} --min-cluster-size={CONFIG["MinClusterSize"]} --n-clusters={CONFIG["N_CLUSTERS"]} --mask={self.spname}/mask-small.png {self.spname}/embeddings-gene.pickle {self.spname}/clusters-gene/')
        run_with_single_thread((f'python {ROOT}/istar/cluster.py ' 
                        f'--filter-size={CONFIG["FilterSize"]} '
                        f'--min-cluster-size={CONFIG["MinClusterSize"]} '
                        f'--n-clusters={CONFIG["N_CLUSTERS"]} '
                        f'--mask={self.spname}/mask-small.png '
                        f'{self.spname}/embeddings-gene.pickle '
                        f'{self.spname}/clusters-gene/'))
        # differential analysis by clusters
        execute_cmd(f'python {ROOT}/istar/aggregate_imputed.py {self.spname}')
        execute_cmd(f'python {ROOT}/istar/reorganize_imputed.py {self.spname}')
        execute_cmd(f'python {ROOT}/istar/differential.py {self.spname}')
        # visualize spot-level gene expression data
        execute_cmd(f'python {ROOT}/istar/plot_spots.py {self.spname}')
        logger.info(f"[{self.spname}] cluster done.")
    

    @timer
    def ln_outs(self) -> None:
        logger.info(f"[{self.spname}] Running ln_outs step.")
        # outs
        mkdir(f'{self.spname}/outs/')
        #os.chdir(f'{self.spname}/outs/') 
        execute_cmd(f'cd {self.spname}/outs/ && ln -s ../clusters-gene/ .')
        execute_cmd(f'cd {self.spname}/outs/ && ln -s ../spots/ .')
        #execute_cmd(f'cd {self.spname}/outs/ && ln -s ../cnts-super-plots/ .')

    @timer
    def run(self) -> None:
        step_order = ['input','preprocess','impute','cluster', 'ln_outs']
        for step in step_order:
            if step in self.step:
                getattr(self, step)()


def parse_mapfile(mapfile, step):
    df_mapfile = pd.read_csv(mapfile, sep=r'\s+')
    df_mapfile['step'] = step
    dir_list = df_mapfile['dir']
    image_list = df_mapfile['image']
    spname_list = df_mapfile['spname']
    swap_pos_list = df_mapfile['swap_pos']
    step_list = df_mapfile['step']
    return dir_list, image_list, spname_list, swap_pos_list, step_list

def run_single(dir, image, spname, swap_pos, step):
    try:
        runner = IStar(dir, image, spname, swap_pos, step)
        runner.run()
        logger.info(f'Completed: {spname}')
    except Exception as e:
        logger.error(f'JOB FAILED: {spname}',
              file=sys.stderr)
        traceback.print_exc()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mapfile', help='tsv:dir<>image<>spname<>swap_pos', required=True)
    parser.add_argument('--step', default='input,preprocess,impute,cluster,ln_outs', help='comma-separated step')
    parser.add_argument('--threads', type=int, default=1, help='thread pool size')
    args = parser.parse_args()

    dir_list, image_list, spname_list, swap_pos_list, step_list = parse_mapfile(args.mapfile, args.step)
    
    logger.info(f"Starting pipeline with {len(dir_list)} samples and {args.threads} thread(s).")
    with ThreadPoolExecutor(max_workers = args.threads) as executor:
        executor.map(run_single, dir_list, image_list, spname_list, swap_pos_list, step_list)


if __name__ == '__main__':
    main()