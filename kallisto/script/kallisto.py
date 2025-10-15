#!/usr/bin/env python

import argparse
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
import glob
import scanpy as sc
import numpy as np

THREADS = 16

class Kallisto():
    def __init__(self, library_id, fq_dir, sample, match_dir, workflow, host_species):
        '''
        '''
        self.library_id = library_id
        self.fq_dir = fq_dir
        self.sample = sample
        self.fq1 = glob.glob(f"{fq_dir}/{library_id}*1*gz")[0]
        self.fq2 = glob.glob(f"{fq_dir}/{library_id}*2*gz")[0]
        self.match_bc = glob.glob(f"{match_dir}/*/*filtered*/barcodes.tsv*")[0]

        if workflow == "sgr":
            self.x = '0,0,9,0,25,34,0,50,59:0,60,72:1,0,0'
            self.w = '/SGRNJ06/randd/USER/zhouyiqi/work/analysis/kb_python/test/1769k-GEXSCOPE-V2.txt'
        
        if host_species == "human":
            self.index = "/SGRNJ06/randd/USER/wangjingshen/rd_project/kallisto/reference/pre_computed/palmdb_human_dlist_cdna_dna.idx"
        if host_species == "mouse":
            self.index = "/SGRNJ06/randd/USER/wangjingshen/rd_project/kallisto/reference/pre_computed/palmdb_mouse_dlist_cdna_dna.idx"


    def run_kb(self):
        cmd = (
            f'kb count '
            f'--verbose '
            f'--aa '
            f'-i {self.index} '
            f'-g /SGRNJ06/randd/USER/wangjingshen/rd_project/kallisto/reference/PalmDB/palmdb_clustered_t2g.txt '
            f'-x {self.x} '
            f'-w {self.w} '
            f'--parity single '
            f'-o {self.sample} '
            f'-t {THREADS} '
            f'--h5ad '
            f'{self.fq1} {self.fq2}'
        )
        subprocess.check_call(cmd, shell=True)

    def h5ad_2_df(self):
        '''
        match bc
        '''
        df = sc.read_h5ad(f"{self.sample}/counts_unfiltered/adata.h5ad")
        columns_with_positive = df.var_names[((df.X > 0).sum(axis=0) > 0).tolist()]     # del undetected virus;  reduce data size
        df = df[:, columns_with_positive].to_df()

        # match bc
        bc = pd.read_csv(self.match_bc, header=None).iloc[:,0]
        bc = bc.apply(lambda x:x.replace("_",""))  # del _
        bc = bc.reset_index()
        bc.columns = ["id","barcode"]
        df_merge = df.merge(bc, how="right", on="barcode")
        df_merge = df_merge.set_index("barcode")
        df_merge = df_merge[df.columns[1:]]
        df_merge = df_merge.fillna(0)     
        df = df_merge.loc[:, df_merge.sum() > 0] # del undetected virus

        # trans id to taxonomy
        ID_to_taxonomy = pd.read_csv("/SGRNJ06/randd/USER/wangjingshen/rd_project/kallisto/reference/PalmDB/ID_to_taxonomy_mapping.csv", index_col=0)
        ID_to_taxonomy = ID_to_taxonomy.loc[ df.columns,:]
        #print(df.columns.equals(ID_to_taxonomy.index))
        df.columns = ID_to_taxonomy.apply(lambda row: '|'.join(row.astype(str)), axis=1)
        df = df.reset_index()
        df.to_csv(f"{self.sample}/counts_unfiltered/virus_df.tsv", sep="\t", index = False)

    def run(self):
        self.run_kb()
        self.h5ad_2_df()

def parse_mapfile(mapfile, workflow, host_species):
    df_mapfile = pd.read_csv(mapfile, sep='\s+', header=None)
    library_id_list = df_mapfile.iloc[:,0]
    fq_dir_list = df_mapfile.iloc[:,1]
    sample_list = df_mapfile.iloc[:,2]
    match_dir_list = df_mapfile.iloc[:,3]

    workflow_list = np.repeat(workflow, df_mapfile.shape[0])
    host_species_list = np.repeat(host_species, df_mapfile.shape[0])
    
    return library_id_list, fq_dir_list, sample_list, match_dir_list, workflow_list, host_species_list


def run_single(library_id, fq_dir, sample, match_dir, workflow, host_species):
    runner = Kallisto(library_id, fq_dir, sample, match_dir, workflow, host_species)
    runner.run()


if __name__ == '__main__':
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--mapfile', help='mapfile', required=True)
    parser.add_argument('--workflow', help='workflow', required=True)
    parser.add_argument('--host_species', required=True)

    # add version
    parser.add_argument('--version', action='version', version='1.0')
    
    args = parser.parse_args()
    library_id_list, fq_dir_list, sample_list, match_dir_list, workflow_list, host_species_list = parse_mapfile(args.mapfile, args.workflow, args.host_species)

    with ProcessPoolExecutor(max_workers = 1) as executor:
        executor.map(run_single, library_id_list, fq_dir_list, sample_list, match_dir_list, workflow_list, host_species_list)