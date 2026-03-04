import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc

import os
import cv2
import json

from pathlib import Path
from PIL import Image



def cut_visium_spots(celescope_space, out_dir='spots_image'):
    """
    Based on the hires graph of adata.uns['spatial'] and the coordinates of obsm['spatial'],
    cut patches with a physical diameter of 55 µm, and save them as png files
    """
    adata = sc.read_visium(celescope_space)
    adata.obsm["spatial"] = adata.obsm["spatial"].astype(float)    # fix 10X str

    img = adata.uns['spatial'][list(adata.uns['spatial'].keys())[0]]['images']['hires']
    img = (img * 255).astype(np.uint8)           # 0-1 → 0-255
    img = cv2.cvtColor(img, cv2.COLOR_RGB2BGR)

    with open(f'{celescope_space}/spatial/scalefactors_json.json') as f:
        sf = json.load(f)
    tissue_hires_scalef = sf['tissue_hires_scalef']

    centers = adata.obsm["spatial"] * tissue_hires_scalef
    spot_df = pd.DataFrame(centers, columns=['x', 'y'], index=adata.obs.index)
    spot_df = spot_df.join(adata.obs[['array_row', 'array_col']]).astype(float)   # fix

    fixed_row = spot_df['array_row'].iloc[0]
    row_df = spot_df[spot_df['array_row'] == fixed_row]
    min_col, max_col = row_df['array_col'].min(), row_df['array_col'].max()
    min_x, max_x = row_df['x'].min(), row_df['x'].max()
    px_per_um = (max_x - min_x) / ((max_col - min_col) / 2) / 100
    radius = int(px_per_um * 55 / 2)

    out_dir = Path(out_dir)
    out_dir.mkdir(exist_ok=True)
    
    invalid_spots = []
    
    h, w = img.shape[:2]
    for idx, row in spot_df.iterrows():
        cx, cy = int(row.x), int(row.y)
        x1 = max(cx - radius, 0)
        y1 = max(cy - radius, 0)
        x2 = min(cx + radius, w)
        y2 = min(cy + radius, h)
        patch = img[y1:y2, x1:x2]

        if patch.size <= 0:
            print(f"Warning: Spot {idx} produces empty patch, skipping")
            print(radius)
            print((h,w))
            print((x1,x2))
            print((y1,y2))
            invalid_spots.append(idx)
            continue

        cv2.imwrite(str(out_dir / f'{idx}_hires.png'), patch)

        # check
        #spot_bc, spot_row, spot_col = idx, int(row.array_row), int(row.array_col)
        #spot_name = f'r{spot_row}_c{spot_col}'
        #cv2.imwrite(str(out_dir / f'{spot_name}.png'), patch)

    print(f'{len(spot_df)} spots saved to {out_dir}')

    return(invalid_spots)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--celescope_space', help='celescope space analysis path', required=True)
    parser.add_argument('--out_dir', default="spots_image", help='outdir')
    args = parser.parse_args()

    cut_visium_spots(args.celescope_space, out_dir = args.out_dir)


if __name__ == '__main__':
    main()