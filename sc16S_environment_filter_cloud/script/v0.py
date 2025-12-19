import pandas as pd
import argparse
import subprocess
from concurrent.futures import ProcessPoolExecutor

N_threads = 4


def environment_filter(fj_dir, environment_genus, sample):
    cmd1 = f'mkdir -p {sample}/outs/'
    cmd2 = f'mkdir -p {sample}/04.analysis_pathseq/'
    subprocess.check_call(cmd1, shell = True)
    subprocess.check_call(cmd2, shell = True)

    raw_matrix_file = f"{fj_dir}/call-pathseq_cellranger3/execution/{sample}_EmptyDrops_CR_raw_UMI_matrix.tsv.gz"
    filter_matrix_file = f"{sample}/outs/{sample}_filter_UMI_matrix.tsv.gz"

    df = pd.read_csv(raw_matrix_file, sep="\t", index_col=0)
    environment_genus = pd.read_csv(environment_genus, header = None).iloc[:0]

    df_filter = df[~df.index.isin(environment_genus)]  # rm environment
    df_filter.to_csv(filter_matrix_file, sep="\t")


def report(sample, match_dir):
    '''
    report
    '''
    cmd = (
        f'celescope pathseq analysis_pathseq '
        f'--outdir {sample}/04.analysis_pathseq '
        f'--sample {sample} '
        f'--thread {N_threads} '
        f'--umi_matrix_file {sample}/outs/{sample}_filter_UMI_matrix.tsv.gz '
        f'--match_dir {match_dir} '
    )
    subprocess.check_call(cmd, shell=True)

def rm_dirs(sample):
    '''
    report
    '''
    cmd = f'rm -rf {sample}/04.analysis_pathseq'
    subprocess.check_call(cmd, shell=True)


def parse_mapfile(mapfile):
    '''
    mapfile has four columns, which are fj_dir, environment_genus, sample, match_dir in order.
    '''
    df_mapfile = pd.read_csv(mapfile, sep='\s+', header = None)
    fj_dir_list = df_mapfile.iloc[:,0]
    environment_genus_list = df_mapfile.iloc[:,1]
    sample_list = df_mapfile.iloc[:,2]
    match_dir_list = df_mapfile.iloc[:,3]
    
    return fj_dir_list, environment_genus_list, sample_list, match_dir_list


def run_single(fj_dir, environment_genus, sample, match_dir):
    environment_filter(fj_dir, environment_genus, sample)
    report(sample, match_dir)
    rm_dirs(sample)

if __name__ == '__main__':
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--mapfile', help="mapfile", required=True)
    args = parser.parse_args()

    fj_dir_list, environment_genus_list, sample_list, match_dir_list = parse_mapfile(args.mapfile)    

    with ProcessPoolExecutor(max_workers = 1) as executor:
        executor.map(run_single, fj_dir_list, environment_genus_list, sample_list, match_dir_list)