import os
import argparse
import subprocess
import pandas as pd


def get_fq(mapfile, x, fq):
    '''
    # prefix   fq_path sample_id   (match_dir)
	CUS241108072    /SGRNJ06/DATA04/limsdownload/24_12/2024_12_22/P24110804/SQ241107XIAZHAI/2024-12-22-10   SQ241107XIAZHAI_Target
    CUS241108072    /SGRNJ06/DATA04/limsdownload/24_12/2024_12_25/P24110804/SQ241107XIAZHAI/2024-12-25-1    SQ241107XIAZHAI_Target
    # x: row of mapfile
    # fq: R1 or R2  
    '''
    fq = f'{mapfile.fq_path.iloc[x]}/{mapfile.prefix.iloc[x]}*{fq}*gz'
    return(fq)

def merge_fq(mapfile, fq, merge_outdir):
    fq_join = ' '.join([get_fq(mapfile, x, fq) for x in range(mapfile.shape[0])])
    outdir = f'{merge_outdir}/{mapfile.sample_id.iloc[0]}_fq/'
    if not os.path.exists(outdir):
        os.system(f"mkdir -p {outdir}")
    cmd = f'cat {fq_join} > {outdir}/{mapfile.prefix.iloc[0]}_{fq}.fastq.gz'
    subprocess.check_call(cmd, shell=True)

def get_fq_path(mapfile, x, merge_outdir):
    merge_path = f'{merge_outdir}/{mapfile.sample_id.iloc[x]}_fq/'
    return(merge_path)

def generate_mapfile(mapfile, merge_outdir, mapfile_outdir):
    mapfile_merge = mapfile.drop_duplicates(subset =['sample_id'])
    mapfile_merge_copy = mapfile_merge.copy()
    mapfile_merge_copy.fq_path = [get_fq_path(mapfile_merge, x, merge_outdir) for x in range(mapfile_merge.shape[0])]
    mapfile_merge_copy.to_csv(f'{mapfile_outdir}/mapfile_merge', sep="\t", index = False, header = False)


def main():
    parsers = argparse.ArgumentParser()
    parsers.add_argument('--mapfile', help='mapfile', required=True)
    parsers.add_argument('--merge_outdir', help='merge outdir', required=True)
    parsers.add_argument('--new_mapfile', help='mapfile outdir', action = "store_true")
    args = parsers.parse_args()

    mapfile = pd.read_csv(args.mapfile, sep="\s+", header=None)
    if(mapfile.shape[1]==3):
        mapfile.columns = ["prefix", "fq_path", "sample_id"]
    if(mapfile.shape[1]==4):
        mapfile.columns = ["prefix", "fq_path", "sample_id", "match_dir"]

    for sample in mapfile.sample_id.unique():
        mapfile_subset = mapfile.loc[mapfile.sample_id==sample]
        merge_fq(mapfile_subset, "R1", args.merge_outdir)  # R1
        merge_fq(mapfile_subset, "R2", args.merge_outdir)  # R2

    if(args.new_mapfile == True):
        generate_mapfile(mapfile, args.merge_outdir, f'{args.merge_outdir}/../')

if __name__ == '__main__':
    main()