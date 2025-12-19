#!/usr/bin/env python

import argparse
import os
import subprocess
import sys
import pandas as pd
from concurrent.futures import ProcessPoolExecutor

from environment_filter import environment_filter

N_threads = 4


class Environment_filter():
    def __init__(self, fj_dir, environment_genus, sample, match_dir, step):
        '''
        '''
        self.fj_dir = fj_dir
        self.environment_genus = environment_genus
        self.sample = sample
        self.match_dir = match_dir
        self.step = args.step

    def ln_mkdir(self):
        cmd1 = f"mkdir -p {self.sample}/"
        cmd2 = f"ln -s {self.fj_dir}/{self.sample}/00.sample/  {self.sample}"
        cmd3 = f"ln -s {self.fj_dir}/{self.sample}/01.starsolo/  {self.sample}"
        cmd4 = f"ln -s {self.fj_dir}/{self.sample}/02.pathseq/  {self.sample}"
        cmd5 = f"ln -s {self.fj_dir}/{self.sample}/03.count_pathseq/  {self.sample}"
        cmd6 = f"ln -s {self.fj_dir}/{self.sample}/04.analysis_pathseq/  {self.sample}"
        cmd7 = f"cp -r {self.fj_dir}/{self.sample}/outs/  {self.sample}"

        subprocess.check_call(cmd1, shell=True)
        subprocess.check_call(cmd2, shell=True)
        subprocess.check_call(cmd3, shell=True)
        subprocess.check_call(cmd4, shell=True)
        subprocess.check_call(cmd5, shell=True)
        subprocess.check_call(cmd6, shell=True)
        subprocess.check_call(cmd7, shell=True)

    def filter(self):
        '''
        '''
        environment_filter(self.environment_genus, self.sample)

    def report(self):
        '''
        report
        '''
        cmd = (
            f'celescope pathseq analysis_pathseq '
            f'--outdir {self.sample}/04.analysis_pathseq '
            f'--sample {self.sample} '
            f'--thread {N_threads} '
            f'--umi_matrix_file {self.sample}/outs/{self.sample}_filter_UMI_matrix.tsv.gz '
            f'--match_dir {self.match_dir} '
        )
        subprocess.check_call(cmd, shell=True)


    def run(self):
        if 'ln_mkdir' in self.step:
            self.ln_mkdir()
        if 'filter' in self.step:
            self.filter()
        if 'report' in self.step:
            self.report()

def parse_mapfile(mapfile, step):
    '''
    mapfile has four columns, which are fj_dir, environment_genus, sample, match_dir in order.
    '''
    df_mapfile = pd.read_csv(mapfile, sep='\s+', header = None)
    df_mapfile.columns = ["fj_dir", "environment_genus", "sample", "match_dir"]
    df_mapfile['step'] = step

    fj_dir_list = df_mapfile['fj_dir']
    environment_genus_list = df_mapfile['environment_genus']
    sample_list = df_mapfile['sample']
    match_dir_list = df_mapfile['match_dir']
    step_list = df_mapfile['step']
    
    return fj_dir_list, environment_genus_list, sample_list, match_dir_list, step_list


def run_single(fj_dir, environment_genus, sample, match_dir, step):
    runner = Environment_filter(fj_dir, environment_genus, sample, match_dir, step)
    runner.run()


if __name__ == '__main__':
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--mapfile', help='mapfile', required=True)
    parser.add_argument('--step')

    # add version
    parser.add_argument('--version', action='version', version='1.0')
    args = parser.parse_args()

    if(args.step is None): # None run all steps
        args.step = 'ln_mkdir,filter,report'

    fj_dir_list, environment_genus_list, sample_list, match_dir_list, step_list = parse_mapfile(args.mapfile, args.step)    

    with ProcessPoolExecutor(max_workers = 1) as executor:
        executor.map(run_single, fj_dir_list, environment_genus_list, sample_list, match_dir_list, step_list)
