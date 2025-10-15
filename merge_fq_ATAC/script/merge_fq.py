import os
import argparse
import subprocess
import pandas as pd


def get_fq(mapfile, x, fq):
    fq = f'{mapfile.fastq_path.iloc[x]}/{mapfile.prefix.iloc[x]}*{fq}*gz'
    return(fq)

def merge_fq(mapfile, fq, merge_outdir):
    outdir = f'{merge_outdir}/{mapfile.outdir.iloc[0]}_fq/'
    if not os.path.exists(outdir):
        os.system(f"mkdir -p {outdir}")

    fq_join = ' '.join([get_fq(mapfile, x, fq) for x in range(mapfile.shape[0])])
    cmd = f'cat {fq_join} > {outdir}/{mapfile.prefix.iloc[0]}_{fq}.fastq.gz'
    print(cmd)
    subprocess.check_call(cmd, shell=True)

def get_fastq_path(mapfile, x, merge_outdir):
    merge_path = f'{merge_outdir}/{mapfile.outdir.iloc[x]}_fq/'
    return(merge_path)

def generate_mapfile(mapfile, merge_outdir, mapfile_outdir):
    mapfile_merge = mapfile.drop_duplicates(subset =['outdir'])
    mapfile_merge_copy = mapfile_merge.copy()
    mapfile_merge_copy.fastq_path = [get_fastq_path(mapfile_merge, x, merge_outdir) for x in range(mapfile_merge.shape[0])]
    mapfile_merge_copy.to_csv(f'{mapfile_outdir}/mapfile_merge', sep="\t", index = False, header = True)


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--mapfile', help='mapfile', required=True)
    parsers.add_argument('--merge_outdir', help='merge outdir', required=True)
    parsers.add_argument('--mapfile_outdir', help='mapfile outdir', required=True)
    args = parsers.parse_args()

    mapfile = pd.read_csv(args.mapfile, sep="\s+", header=0)
    #if(mapfile.shape[1]==3):
    #    mapfile.columns = ["prefix", "fastq_path", "outdir"]
    #if(mapfile.shape[1]==4):
    #    mapfile.columns = ["prefix", "fastq_path", "outdir", "species"]
    print(mapfile.columns)  

    merge_outdir = args.merge_outdir
    mapfile_outdir = args.mapfile_outdir

    for sample in mapfile.outdir.unique():
        mapfile_subset = mapfile.loc[mapfile.outdir==sample]
        #print(mapfile_subset)
        merge_fq(mapfile_subset, "R1", merge_outdir)  # R1
        merge_fq(mapfile_subset, "R2", merge_outdir)  # R2

    generate_mapfile(mapfile, merge_outdir, mapfile_outdir)

if __name__ == '__main__':
    main()