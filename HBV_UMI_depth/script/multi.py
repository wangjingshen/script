import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
from get_virus_bam import Bam_virus

class multi():
    def __init__(self, fj_path, otsu_min_support_read, outdir, name):
        self.fj_path = fj_path
        self.otsu_min_support_read = otsu_min_support_read
        self.outdir = outdir
        self.name = name

    def run(self):
        Bam_virus(fj_path, otsu_min_support_read, outdir, name)

def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, sep='\t')
    fj_path_list = df_mapfile['fj_path']
    otsu_min_support_read_list = df_mapfile['otsu_min_support_read']
    outdir_list = df_mapfile['outdir']
    name_list = df_mapfile['name']
    
    return fj_path_list, otsu_min_support_read_list, outdir_list, name_list

def run_single(fj_path, otsu_min_support_read, outdir, name):
    runner = multi(fj_path, otsu_min_support_read, outdir, name)
    runner.run()

def main():
    mapfile = sys.argv[1]
    fj_path_list, otsu_min_support_read_list, outdir_list, name_list = parse_mapfile(mapfile)
    print(fj_path_list)
    with ProcessPoolExecutor(max_workers = 1) as executor:
        executor.map(run_single, fj_path_list, otsu_min_support_read_list, outdir_list, name_list)


if __name__ == '__main__':
    main()
