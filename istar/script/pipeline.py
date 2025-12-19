import sys
import pandas as pd
import argparse
import os
import glob
import subprocess
from concurrent.futures import ProcessPoolExecutor

ROOT = os.path.dirname(os.path.abspath(__file__))

DEVICE = "cpu"  # "cuda" or "cpu"

RADIUS_RAW = 20
PixelSizeRaw = 0.25
PixelSize = 0.5
N_HVG = 1000
EPOCHS = 400
FilterSize = 8 
MinClusterSize = 20
N_CLUSTERS = 10

class IStar:
    def __init__(self, dir, image, pos, spname, step):
        self.mtx = glob.glob(f"{dir}/outs/filtered")[0]
        self.image = image #glob.glob(f'{dir}/spatial/figure/*jpg')[0]
        self.pos = pos
        self.spname = spname+"/"  #istar prefix
        self.step = step.strip().split(',')
    
    def input(self):
        '''
        generate_input
        '''
        if not os.path.exists(self.spname):
            os.system(f"mkdir -p {self.spname}")
        
        cmd1 = (f'/SGRNJ/Public/Software/conda_env/r4.1_env/bin/Rscript {ROOT}/gen_cnts_locs.R '
                f'--mtx {self.mtx} '
                f'--pos {self.pos} '
                f'--spname {self.spname} '
                )
        cmd2 = f'cp -r {ROOT}/../data/checkpoints/ .'
        cmd3 = f'echo {RADIUS_RAW} > {self.spname}/radius-raw.txt'
        cmd4 = f'echo {PixelSizeRaw} > {self.spname}/pixel-size-raw.txt'
        cmd5 = f'echo {PixelSize} > {self.spname}/pixel-size.txt'
        cmd6 = f'cp {self.image} {self.spname}/he-raw.png'

        subprocess.check_call(cmd1, shell = True)
        subprocess.check_call(cmd2, shell = True)
        subprocess.check_call(cmd3, shell = True)
        subprocess.check_call(cmd4, shell = True)
        subprocess.check_call(cmd5, shell = True)
        subprocess.check_call(cmd6, shell = True)


    def preprocess(self):
        '''
        preprocess_image
        '''
        cmd1 = f'python {ROOT}/istar/rescale.py {self.spname}/ --image'
        cmd2 = f'python {ROOT}/istar/preprocess.py {self.spname} --image'
        # extract histology features
        cmd3 = f'python {ROOT}/istar/extract_features.py {self.spname} --device={DEVICE}'
        # auto detect tissue mask
        cmd4 = f'python {ROOT}/istar/get_mask.py {self.spname}/embeddings-hist.pickle {self.spname}/mask-small.png'
        # select most highly variable genes to predict
        cmd5 = f'python {ROOT}/istar/select_genes.py --n-top={N_HVG} {self.spname}/cnts.tsv {self.spname}/gene-names.txt'
        # rescale coordinates and spot radius
        cmd6 = f'python {ROOT}/istar/rescale.py {self.spname} --locs --radius'

        subprocess.check_call(cmd1, shell = True)
        subprocess.check_call(cmd2, shell = True)
        subprocess.check_call(cmd3, shell = True)
        subprocess.check_call(cmd4, shell = True)
        subprocess.check_call(cmd5, shell = True)
        subprocess.check_call(cmd6, shell = True)


    def impute(self):
        '''
        predict super-resolution gene expression
        '''
        # train gene expression prediction model and predict at super-resolution
        cmd1 = f'python {ROOT}/istar/impute.py {self.spname} --epochs={EPOCHS} --device={DEVICE}' 
        # visualize imputed gene expression
        cmd2 = f'python {ROOT}/istar/plot_imputed.py {self.spname}'

        subprocess.check_call(cmd1, shell = True)
        subprocess.check_call(cmd2, shell = True)


    def cluster(self):
        # segment image by gene features
        cmd1 = f'python {ROOT}/istar/cluster.py --filter-size={FilterSize} --min-cluster-size={MinClusterSize} --n-clusters={N_CLUSTERS} --mask={self.spname}/mask-small.png {self.spname}/embeddings-gene.pickle {self.spname}/clusters-gene/'
        # differential analysis by clusters
        cmd2 = f'python {ROOT}/istar/aggregate_imputed.py {self.spname}'
        cmd3 = f'python {ROOT}/istar/reorganize_imputed.py {self.spname}'
        cmd4 = f'python {ROOT}/istar/differential.py {self.spname}'
        # visualize spot-level gene expression data
        cmd5 = f'python {ROOT}/istar/plot_spots.py {self.spname}'

        subprocess.check_call(cmd1, shell = True)
        subprocess.check_call(cmd2, shell = True)
        subprocess.check_call(cmd3, shell = True)
        subprocess.check_call(cmd4, shell = True)
        subprocess.check_call(cmd5, shell = True)


    def run(self):
        if 'input' in self.step:
            self.input()
        if 'preprocess' in self.step:
            self.preprocess()
        if 'impute' in self.step:
            self.impute()
        if 'cluster' in self.step:
            self.cluster()


def parse_mapfile(mapfile, step):
    df_mapfile = pd.read_csv(mapfile, sep='\s+')
    if(step is None):
        step = 'input,preprocess,impute,cluster'
    df_mapfile['step'] = step
    
    dir_list = df_mapfile['dir']
    image_list = df_mapfile['image']
    pos_list = df_mapfile['pos']
    spname_list = df_mapfile['spname']
    step_list = df_mapfile['step']

    return dir_list, image_list, pos_list, spname_list, step_list

def run_single(dir, image, pos, spname, step):
    runner = IStar(dir, image, pos, spname, step)
    print(runner)
    runner.run()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mapfile', help='mapfile', required=True)
    parser.add_argument('--step', help='step')
    parser.add_argument('--version', action='version', version='1.0')
    args = parser.parse_args()

    dir_list, image_list, pos_list, spname_list, step_list = parse_mapfile(args.mapfile, args.step)
    print(step_list)
    with ProcessPoolExecutor(max_workers = 1) as executor:
        executor.map(run_single, dir_list, image_list, pos_list, spname_list, step_list)


if __name__ == '__main__':
    main()