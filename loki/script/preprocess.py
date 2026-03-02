import sys
import argparse
from pathlib import Path

from get_coord import get_coord
from h5toh5ad import h5toh5ad
from cut_visium_spots import cut_visium_spots
from encode import encode_space_text, encode_space_image, encode_sc_text, make_finetune

dev_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(dev_root))
from utils.utils import mkdir, logger, execute_cmd, timer

class Loki_preprocess:
    def __init__(self, dir: Path, spname: str, hk_genes: str, sc_h5ad: str, step: str):
        self.dir = dir
        self.spname = spname
        self.space_input = f'{self.spname}/space_input/'
        self.hk_genes = hk_genes
        self.sc_h5ad = sc_h5ad
        self.loki_input = f'{self.spname}/loki_input/'
        self.step = step.strip().split(',')

    def space_preprocess(self) -> None:
        '''
        generate_input
        '''
        logger.info(f"[{self.spname}] Running space_preprocess step.")
        mkdir(self.space_input)
        execute_cmd((f'cp {self.dir}/outs/filtered_feature_bc_matrix.h5 {self.space_input}'))
        execute_cmd((f'cp -r {self.dir}/outs/spatial {self.space_input}'))
        execute_cmd((f'mv {self.space_input}/spatial/positions_list.csv {self.space_input}/spatial/tissue_positions_list.csv'))
        logger.info(f"[{self.spname}] space_preprocess done.")

    def h5toh5ad(self) -> None:
        '''
        h5 to h5ad
        '''
        h5toh5ad(self.space_input)
    
    def loki_preprocess(self) -> None:
        '''
        preprocess_image
        '''
        logger.info(f"[{self.spname}] Running loki_preprocess step.")
        mkdir(self.loki_input)
        if self.sc_h5ad is not None:
            encode_sc_text(self.sc_h5ad, self.hk_genes, self.loki_input)
        get_coord(self.space_input)
        cut_visium_spots(self.space_input, f'{self.loki_input}/spots_images/')
        encode_space_text(self.space_input, self.hk_genes, self.loki_input)
        encode_space_image(self.space_input, f'{self.space_input}/image_coord.csv', self.loki_input)

        logger.info(f"[{self.spname}] loki_preprocess done.")
    
    def loki_preprocess(self) -> None:
        '''
        preprocess_image
        '''
        logger.info(f"[{self.spname}] Running loki_preprocess step.")
        mkdir(self.loki_input)
        if self.sc_h5ad is not None:
            encode_sc_text(self.sc_h5ad, self.hk_genes, self.loki_input)
        get_coord(self.space_input)
        cut_visium_spots(self.space_input, f'{self.loki_input}/spots_images/')
        encode_space_text(self.space_input, self.hk_genes, self.loki_input)
        encode_space_image(self.space_input, f'{self.space_input}/image_coord.csv', self.loki_input)

        logger.info(f"[{self.spname}] loki_preprocess done.")

    @timer
    def run(self) -> None:
        step_order = ['space_preprocess','h5toh5ad','loki_preprocess']
        for step in step_order:
            if step in self.step:
                getattr(self, step)()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', help='celescope sapce dir', required=True)
    parser.add_argument('--spname', help='spname', required=True)
    parser.add_argument('--hk_genes', help='housekeeping_genes', required=True)
    parser.add_argument('--sc_h5ad', help='sc h5ad', default=None)
    parser.add_argument('--step', default='space_preprocess,h5toh5ad,loki_preprocess', help='comma-separated step')
    args = parser.parse_args()

    runner = Loki_preprocess(args.dir, args.spname, args.hk_genes, args.sc_h5ad, args.step)
    runner.run()


if __name__ == '__main__':
    main()