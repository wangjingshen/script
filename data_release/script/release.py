import os
import csv
import sys
import argparse
import json
import multiprocessing
import pathlib

from tqdm import tqdm
from collections import defaultdict


def read_json(data_json):
    with open(data_json, 'r') as fh:
        data = json.load(fh)
    return data

def cp_fastq(params):
    sample,prefix_dict,library_dict,outdir = params
    os.system(f'mkdir -p {outdir}/{sample}')
    folders = library_dict[sample]
    for folder in folders:
        folder = os.path.abspath(folder)
        for root, _, files in os.walk(folder):
            for _file in files:
                if not _file.endswith(".gz"):
                    continue
                for prefix in set(prefix_dict[sample]):
                    if prefix+"_" in _file:
                        path = os.path.join(root, _file)
                        cmd = f'cp {path} {outdir}/{sample}'
                        os.system(cmd)

def cp_fastqs(mapfile_file,outdir):
    prefix_dict = defaultdict(list)
    library_dict = defaultdict(list)
    with open(mapfile_file, 'r') as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\s+')
        for row in reader:
            prefix_dict[row['sample']].append(row['prefix'])
            library_dict[row['sample']].append(row['library'])
    param_list = []
    for sample in prefix_dict:
        param_list.append((sample,prefix_dict,library_dict,outdir))
    with multiprocessing.Pool(len(prefix_dict)) as p:
        list(tqdm(p.imap(cp_fastq,param_list),total=len(param_list),unit_scale = True,ncols = 70,file = sys.stdout,desc='Cp fastq '))
    p.close()
    p.join()

def cp_matrixs(matrix_type,in_outdir,outdir):
    if matrix_type == 'rna':
        for sample in os.listdir(in_outdir):
            if pathlib.Path(f'{in_outdir}/{sample}/outs').is_dir():
                os.system(f'mkdir -p {outdir}/matrix/rna/{sample}')
                matrix_path_f = f'{in_outdir}/{sample}/outs/filtered'
                #matrix_path_r = f'{in_outdir}/{sample}/outs/raw'
                cmd = f'cp -r {matrix_path_f} {outdir}/matrix/rna/{sample}'
                os.system(cmd)
                #cmd = f'cp -r {matrix_path_r} {outdir}/matrix/rna/{sample}'
                #os.system(cmd)
    elif matrix_type == 'tag':
        for sample in os.listdir(in_outdir):
            if pathlib.Path(f'{in_outdir}/{sample}/03.count_tag').is_dir():
                os.system(f'mkdir -p {outdir}/matrix/tag/{sample}')
                matrix_path = f'{in_outdir}/{sample}/03.count_tag/{sample}_umi_tag.tsv'
                cmd = f'cp {matrix_path} {outdir}/matrix/tag/{sample}'
                os.system(cmd)
    elif matrix_type == 'ebv':
        for sample in os.listdir(in_outdir):
            if pathlib.Path(f'{in_outdir}/{sample}/09.count').is_dir():
                os.system(f'mkdir -p {outdir}/matrix/ebv/{sample}')
                matrix_path = f'{in_outdir}/{sample}/09.count/{sample}_virus_matrix'
                cmd = f'cp -r {matrix_path} {outdir}/matrix/ebv/{sample}'
                os.system(cmd)
    elif matrix_type == 'tcr':
        for sample in os.listdir(in_outdir):
            if pathlib.Path(f'{in_outdir}/{sample}/05.match').is_dir():
                os.system(f'mkdir -p {outdir}/matrix/tcr/{sample}')
                matrix_path = f'{in_outdir}/{sample}/05.match/matched_productive_contig_annotations.csv'
                cmd = f'cp {matrix_path} {outdir}/matrix/tcr/{sample}'
                os.system(cmd)
    elif matrix_type == 'bcr':
        for sample in os.listdir(in_outdir):
            if pathlib.Path(f'{in_outdir}/{sample}/05.match').is_dir():
                os.system(f'mkdir -p {outdir}/matrix/bcr/{sample}')
                matrix_path = f'{in_outdir}/{sample}/05.match/matched_productive_contig_annotations.csv'
                cmd = f'cp {matrix_path} {outdir}/matrix/bcr/{sample}'
                os.system(cmd)
    elif matrix_type == 'snp':
        for sample in os.listdir(in_outdir):
            if pathlib.Path(f'{in_outdir}/{sample}/09.analysis_snp').is_dir():
                os.system(f'mkdir -p {outdir}/matrix/snp/{sample}')
                matrix_path = f'{in_outdir}/{sample}/09.analysis_snp/{sample}_gt.csv'
                cmd = f'cp {matrix_path} {outdir}/matrix/snp/{sample}'
                os.system(cmd)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='data release. \n')
    parser.add_argument('--mapfile_json',help='input path name file',required=True)
    parser.add_argument('--outdir',help='results dir')
    parser.add_argument('--data_type',help='Types of data to be uploaded.Choose from fastq,matrix and BAM, multiple selections are allowed, separated by ",".')
    parser.add_argument('--rna_outdir',help='rna outdir')
    parser.add_argument('--tag_outdir',help='tag outdir')
    parser.add_argument('--ebv_outdir',help='ebv outdir')
    parser.add_argument('--bcr_outdir',help='bcr outdir')
    parser.add_argument('--tcr_outdir',help='tcr outdir')
    parser.add_argument('--snp_outdir',help='snp outdir')
    
    args = parser.parse_args()

    mapfile_json = args.mapfile_json
    outdir = args.outdir
    data_type = set(args.data_type.split(','))

    data_dict = read_json(mapfile_json)
    if 'fastq' in data_type:
        for tag in data_dict:
            cp_fastqs(data_dict[tag],f'{outdir}/fastq/{tag}')
    if 'rna_matrix' in data_type:
        rna_outdir = args.rna_outdir.split(",")
        print(rna_outdir)
        res = [cp_matrixs('rna', x, outdir) for x in rna_outdir]
    if 'tag_matrix' in data_type:
        tag_outdir = args.tag_outdir.split(",")
        res = [cp_matrixs('tag', x, outdir) for x in tag_outdir]
    if 'ebv_matrix' in data_type:
        ebv_outdir = args.ebv_outdir.split(",")
        res = [cp_matrixs('ebv', x, outdir) for x in ebv_outdir]
    if 'bcr_matrix' in data_type:
        bcr_outdir = args.bcr_outdir.split(",")
        res = [cp_matrixs('bcr', x, outdir) for x in bcr_outdir]
    if 'tcr_matrix' in data_type:
        tcr_outdir = args.tcr_outdir.split(",")
        res = [cp_matrixs('tcr', x, outdir) for x in tcr_outdir]
    if 'snp_matrix' in data_type:
        snp_outdir = args.snp_outdir.split(",")
        res = [cp_matrixs('snp', x, outdir) for x in snp_outdir]

    cmd_tar = f'tar -zcvf {outdir}.tar.gz {outdir}'
    os.system(cmd_tar)
        