import scanpy as sc
import pandas as pd
import numpy as np
from PIL import Image
import json


def h5toh5ad(space_input):
    adata = sc.read_visium(space_input)
    adata.var_names_make_unique()
    adata.obsm["spatial"] = adata.obsm["spatial"].astype(float)    # fix 10X str
    

    adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(adata, target_sum=1e4, inplace=True)
    sc.pp.log1p(adata)
    
    adata.layers['normalised'] = adata.X

    sc.pp.highly_variable_genes(
        adata, 
        layer=None, 
        n_top_genes=None, 
        min_disp=0.5, 
        max_disp=np.inf, 
        min_mean=0.0125, 
        max_mean=3, 
        span=0.3,
        n_bins=20, 
        flavor="seurat", 
        subset=False, 
        inplace=True, 
        batch_key=None, 
        check_values=True
    )


    sc.pp.scale(
        adata, 
        zero_center=True, 
        max_value=10, 
        copy=False, 
        layer=None, 
        obsm=None 
    )

    adata.layers['scaled'] = adata.X

    sc.pp.pca(
        adata,
        n_comps=50,
        zero_center=True,
        svd_solver="auto",
        random_state=0,
        return_info=False,
        use_highly_variable=True,
        dtype="float32",
        copy=False,
        chunked=False,
        chunk_size=None,
    )

    
    sc.pp.neighbors(
        adata,
        n_neighbors=15,
        n_pcs=25,
        use_rep=None,
        knn=True,
        random_state=0,
        method="umap",
        metric="euclidean",
        key_added=None,
        copy=False,
    )

    
    sc.tl.tsne(
        adata,
        n_pcs=25,
        copy=False,
    )

    
    sc.tl.umap(
        adata,
        min_dist=0.5,
        spread=1.0,
        n_components=2,
        maxiter=None,
        alpha=1.0,
        gamma=1.0,
        negative_sample_rate=5,
        init_pos="spectral",
        random_state=0,
        a=None,
        b=None,
        copy=False,
        method="umap",
        neighbors_key=None,
    )


    sc.tl.leiden(
        adata,
        resolution = 0.8,
        key_added="cluster",
    )

    adata.X = adata.layers['normalised']
    adata.write(f'{space_input}/filtered_feature_bc_matrix.h5ad')