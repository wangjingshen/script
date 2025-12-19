import pandas as pd
import argparse
import subprocess
from concurrent.futures import ProcessPoolExecutor
import json

N_threads = 4
Environment_genus = "/SGRNJ06/randd/USER/wangjingshen/script/sc16S_environment_filter_cloud/data/environment_list.tsv"
Data_json = "/SGRNJ06/randd/USER/wangjingshen/script/sc16S_environment_filter_cloud/data/data.json" # data_json_template


def environment_filter(fj_dir, sample):
    cmd1 = f'mkdir -p {sample}/outs/'
    cmd2 = f'mkdir -p {sample}/04.analysis_pathseq/'
    subprocess.check_call(cmd1, shell = True)
    subprocess.check_call(cmd2, shell = True)

    raw_matrix_file = f"{fj_dir}/call-pathseq_cellranger3/execution/{sample}_EmptyDrops_CR_raw_UMI_matrix.tsv.gz"
    filter_matrix_file = f"{sample}/outs/{sample}_filter_UMI_matrix.tsv.gz"

    df = pd.read_csv(raw_matrix_file, sep="\t", index_col=0)
    environment_genus = pd.read_csv(Environment_genus, header = None).iloc[:,0]
    df_filter = df[~df.index.isin(environment_genus)]  # rm environment
    df_filter.to_csv(filter_matrix_file, sep="\t")


def gen_data_json(fj_dir, sample):
    with open(Data_json, 'r', encoding='utf-8') as file:
        data = json.load(file)
    with open(f'{fj_dir}/call-pathseq_cellranger3/execution/.EmptyDrops_CR_metrics.json', 'r', encoding='utf-8') as file:
        metrics = json.load(file)

    for item in data.keys():
        for i in range(len(data[item]['metric_list'])):
            update_var = data[item]['metric_list'][i]['name']
            if(update_var !='Chemistry'):
                update_values = metrics[item][update_var]
                if type(update_values) == str:
                    data[item]['metric_list'][i]['value'] = update_values
                    data[item]['metric_list'][i]['display'] = update_values
                if type(update_values) == int:
                    data[item]['metric_list'][i]['value'] = update_values
                    data[item]['metric_list'][i]['display'] = f"{update_values:,}"
                if type(update_values) == float:
                    data[item]['metric_list'][i]['value'] = update_values
                    data[item]['metric_list'][i]['display'] = str(update_values)+"%"

    with open(f'{sample}/.data.json', 'w', encoding='utf-8') as file:
        json.dump(data, file, ensure_ascii=False, indent=4)


def report(fj_dir, sample, match_dir):
    '''
    report
    '''
    #cmd1 = f'cp {fj_dir}/call-pathseq_cellranger3/execution/{sample}_EmptyDrops_CR_report.html {sample}/{sample}_report.html'
    cmd2 = (
        f'celescope pathseq analysis_pathseq '
        f'--outdir {sample}/04.analysis_pathseq '
        f'--sample {sample} '
        f'--thread {N_threads} '
        f'--umi_matrix_file {sample}/outs/{sample}_filter_UMI_matrix.tsv.gz '
        f'--match_dir {match_dir} '
    )
    #subprocess.check_call(cmd1, shell=True)
    subprocess.check_call(cmd2, shell=True)


def rm_dirs(sample):
    cmd = f'rm -rf {sample}/04.analysis_pathseq'
    subprocess.check_call(cmd, shell=True)


def parse_mapfile(mapfile):
    '''
    mapfile has 3 columns, which are fj_dir, sample, match_dir in order.
    '''
    df_mapfile = pd.read_csv(mapfile, sep='\s+', header = None)
    fj_dir_list = df_mapfile.iloc[:,0]
    sample_list = df_mapfile.iloc[:,1]
    match_dir_list = df_mapfile.iloc[:,2]
    
    return fj_dir_list, sample_list, match_dir_list


def run_single(fj_dir, sample, match_dir):
    environment_filter(fj_dir, sample)
    gen_data_json(fj_dir, sample)
    report(fj_dir, sample, match_dir)
    rm_dirs(sample)

if __name__ == '__main__':
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--mapfile', help="mapfile", required=True)
    args = parser.parse_args()

    fj_dir_list, sample_list, match_dir_list = parse_mapfile(args.mapfile)    

    with ProcessPoolExecutor(max_workers = 1) as executor:
        executor.map(run_single, fj_dir_list, sample_list, match_dir_list)