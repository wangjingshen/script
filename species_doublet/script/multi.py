import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
import pandas as pd

class Multi():
    def __init__(self, data, percent_threshold, mode, name):
        self.data = data
        self.percent_threshold = percent_threshold
        self.mode = mode
        self.name = name
    
    def run(self):
        cmd = (
            f'Rscript /SGRNJ06/randd/USER/wangjingshen/script_dev/species_doublet/script/analysis.R '
            f'--data {self.data} '
            f'--percent_threshold {self.percent_threshold} '
            f'--mode {self.mode} '
            f'--name {self.name} '
        )
        subprocess.check_call(cmd, shell=True)


def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, header= 0, sep="\s+")
    data_list = df_mapfile['data']
    percent_threshold_list = df_mapfile['percent_threshold']
    mode_list = df_mapfile['mode']
    name_list = df_mapfile['name']
    return data_list, percent_threshold_list, mode_list, name_list

def run_single(data, percent_threshold, mode, name):
    runner = Multi(data, percent_threshold, mode, name)
    runner.run()

def main():
    mapfile = sys.argv[1]
    data_list, percent_threshold_list, mode_list, name_list = parse_mapfile(mapfile)

    with ProcessPoolExecutor(max_workers=1) as executor:
        executor.map(run_single, data_list, percent_threshold_list, mode_list, name_list)


if __name__ == '__main__':
    main()