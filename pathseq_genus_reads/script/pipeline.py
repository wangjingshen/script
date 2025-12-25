#!/usr/bin/env python

import argparse
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
import glob
from pathlib import Path
import logging

dev_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(dev_root))
#print(sys.path)
from utils.utils import find_file, mkdir, execute_cmd


logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)


class Pathseq_genus_reads():
    def __init__(self, dir, sample, tax_name, match_dir, outdir):
        self.sample = sample
        self.pathseq_score = find_file(f"{dir}/{sample}/02.pathseq/{sample}_pathseq_score.txt")
        self.pathseq_bam = find_file(f"{dir}/{sample}/02.pathseq/{sample}_pathseq.bam")
        self.downsample_bam = find_file(f"{dir}/{sample}/02.pathseq/{sample}_downsample.bam")
        self.umi_mtx = find_file(f"{dir}/{sample}/outs/{sample}_EmptyDrops_CR_raw_UMI_matrix.tsv.gz")
        self.tax_name = tax_name
        self.match_dir = match_dir
        self.outdir = outdir

        mkdir(self.outdir)


    def get_tax_id(self):
        logger.info(f"[{self.sample}] Running get_tax_id step.")
        execute_cmd(
            f'cut -f 1,2 {self.pathseq_score} | grep {self.tax_name} > {self.outdir}/{self.sample}_tax_id.tsv '
        )
        logger.info(f"[{self.sample}] Running get_tax_id step.")


    def awk(self):
        logger.info(f"[{self.sample}] Running awk step.")
        execute_cmd(
            f'/SGRNJ/Public/Software/conda_env/celescope3.0.0/bin/samtools view {self.downsample_bam} | '
            '''awk -v FS="\t" -v OFS="\t" '{print $1,$17,$18}' >  '''
            f'{self.outdir}/{self.sample}_downsample.tsv'
        )
        execute_cmd(
            f'/SGRNJ/Public/Software/conda_env/celescope3.0.0/bin/samtools view {self.pathseq_bam} | '
            '''awk -v FS="\t" -v OFS="\t" '{print $1,$16}' >  '''
            f'{self.outdir}/{self.sample}_pathseq.tsv'
        )
        logger.info(f"[{self.sample}] Running awk step.")

    def get_tax_reads(self):
        logger.info(f"[{self.sample}] Running get_tax_reads step.")
        execute_cmd(
            f'Rscript /SGRNJ06/randd/USER/wangjingshen/script_dev/pathseq_genus_reads/script/get_tax_reads.R '
            f'--sample {self.sample} '
            f'--pathseq_tsv {self.outdir}/{self.sample}_pathseq.tsv '
            f'--downsample_tsv {self.outdir}/{self.sample}_downsample.tsv '
            f'--umi_mtx {self.umi_mtx} '
            f'--tax_id {self.outdir}/{self.sample}_tax_id.tsv '
            f'--tax_name {self.tax_name} '
            f'--match_dir {self.match_dir} '
            f'--outdir {self.outdir} '
        )
        execute_cmd(f'rm {self.outdir}/{self.sample}_downsample.tsv {self.outdir}/{self.sample}_pathseq.tsv {self.outdir}/{self.name}_tax_id.tsv')
        logger.info(f"[{self.sample}] Running get_tax_reads step.")


    def run(self):
        self.get_tax_id()
        self.awk()
        self.get_tax_reads()


def parse_mapfile(mapfile, tax_name):
    df_mapfile = pd.read_csv(mapfile, sep='\s+')
    df_mapfile['tax_name'] = tax_name

    dir_list = df_mapfile['dir']
    sample_list = df_mapfile['sample']
    tax_name_list = df_mapfile['tax_name']
    match_dir_list = df_mapfile['match_dir']
    outdir_list = df_mapfile['outdir']
    return dir_list, sample_list, tax_name_list, match_dir_list, outdir_list


def run_single(dir, sample, tax_name, match_dir, outdir):
    runner = Pathseq_genus_reads(dir, sample, tax_name, match_dir, outdir)
    runner.run()
    print(f"{sample} done.")


if __name__ == '__main__':
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--mapfile', help='tsv:dir<>image<>pos<>sample', required=True)
    parser.add_argument('--tax_name', help='tax_name', required=True)
    parser.add_argument('--threads', type=int, default=4, help='tax_name')

    args = parser.parse_args()
    dir_list, sample_list, tax_name_list, match_dir_list, outdir_list = parse_mapfile(args.mapfile, args.tax_name)
    with ThreadPoolExecutor(max_workers = args.threads) as executor:
        executor.map(run_single, dir_list, sample_list, tax_name_list, match_dir_list, outdir_list)