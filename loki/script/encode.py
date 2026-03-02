import os
import pandas as pd
import numpy as np
import scanpy as sc
import torch
import loki.utils
import loki.preprocess
from PIL import Image
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
model_path = f'{ROOT}/data/checkpoint.pt'
device = 'cpu'


def encode_space_text(space_input, housekeeping_genes, outdir) -> None:
    '''
    for align, anno, decompose
    '''
    ad = sc.read_h5ad(f'{space_input}/filtered_feature_bc_matrix.h5ad')   # ad = sc.read_10x_h5(f'{space_input}/filtered_feature_bc_matrix.h5')  ## result same
    ad.var_names_make_unique()
    house_keeping_genes = pd.read_csv(housekeeping_genes, index_col = 0)
    top_k_genes_str = loki.preprocess.generate_gene_df(ad, house_keeping_genes)

    model, preprocess, tokenizer = loki.utils.load_model(model_path, device)
    text_embeddings = loki.utils.encode_text_df(model, tokenizer, top_k_genes_str, 'label', device)
    df = pd.DataFrame(text_embeddings.cpu().numpy().T)
    df.to_csv(f'{outdir}/space_text_embeddings.csv', header=False, float_format='%.6f')  # for align, anno
    #np.save(f'{outdir}/text_embeddings.npy', df) 
    torch.save(text_embeddings, f'{outdir}/space_text_embeddings.pt')
    df.columns = ad.obs.index
    df.to_csv(f'{outdir}/space_text_embeddings_header.csv', float_format='%.6f')  # for decompose


def encode_space_image(space_input, coord, outdir) -> None:
    '''
    for align, anno, decompose
    '''
    ad = sc.read_10x_h5(f'{space_input}/filtered_feature_bc_matrix.h5')
    ad.var_names_make_unique()
    coord = pd.read_csv(coord, index_col=0)
    img = Image.open(f'{space_input}/spatial/tissue_hires_image.png')
    img_array = np.asarray(img)

    # segment
    patch_dir = f'{outdir}/spots_images'
    img_list = os.listdir(patch_dir)
    patch_paths = [os.path.join(patch_dir, fn) for fn in img_list]

    model, preprocess, tokenizer = loki.utils.load_model(model_path, device)
    image_embeddings = loki.utils.encode_images(model, preprocess, patch_paths, device)
    df = pd.DataFrame(image_embeddings.cpu().numpy().T)
    df.to_csv(f'{outdir}/space_image_embeddings.csv', header=False, float_format='%.6f')
    torch.save(image_embeddings, f'{outdir}/space_image_embeddings.pt')
    df.columns = img_list
    df.to_csv(f'{outdir}/space_image_embeddings_header.csv', float_format='%.6f')

    ## loki segment
    patch_dir = f'{outdir}/spots_images_loki'
    loki.preprocess.segment_patches(img_array, coord, patch_dir)
    img_list = os.listdir(patch_dir)
    patch_paths = [os.path.join(patch_dir, fn) for fn in img_list]

    model, preprocess, tokenizer = loki.utils.load_model(model_path, device)
    image_embeddings = loki.utils.encode_images(model, preprocess, patch_paths, device)
    df = pd.DataFrame(image_embeddings.cpu().numpy().T)
    df.to_csv(f'{outdir}/space_image_embeddings_loki.csv', header=False, float_format='%.6f')
    torch.save(image_embeddings, f'{outdir}/space_image_embeddings_loki.pt')
    df.columns = img_list
    df.to_csv(f'{outdir}/space_image_embeddings_header_loki.csv', float_format='%.6f')


def encode_sc_text(sc_h5ad, housekeeping_genes, outdir) -> None:
    '''
    for decompose
    '''
    ad = sc.read_h5ad(sc_h5ad)
    ad.var_names_make_unique()
    house_keeping_genes = pd.read_csv(housekeeping_genes, index_col = 0)
    top_k_genes_str = loki.preprocess.generate_gene_df(ad, house_keeping_genes)

    model, preprocess, tokenizer = loki.utils.load_model(model_path, device)
    text_embeddings = loki.utils.encode_text_df(model, tokenizer, top_k_genes_str, 'label', device)
    df = pd.DataFrame(text_embeddings.cpu().numpy().T)
    df.columns = ad.obs.index
    df.to_csv(f'{outdir}/sc_text_embeddings_header.csv', float_format='%.6f')

