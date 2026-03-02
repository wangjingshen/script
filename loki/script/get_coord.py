import os
import pandas as pd
import numpy as np
import scanpy as sc
import json


def get_coord(space_input):
    # coord
    adata = sc.read_visium(space_input)
    adata.obsm["spatial"] = adata.obsm["spatial"].astype(float)    # fix 10X str
    
    with open(f'{space_input}/spatial/scalefactors_json.json') as f:
        sf = json.load(f)
    tissue_hires_scalef = sf['tissue_hires_scalef']

    centers = adata.obsm["spatial"] * tissue_hires_scalef
    spot_df = pd.DataFrame(centers, columns=['x', 'y'], index=adata.obs.index)
    spot_df.columns = ['pixel_x', 'pixel_y']
    #spot_df = spot_df.join(adata.obs[['array_row', 'array_col']])
    spot_df.to_csv(f'{space_input}/image_coord.csv')


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--space_input', help='space input', required=True)

    args = parsers.parse_args()
    get_coord(args.space_dir)

if __name__ == '__main__':
    main()