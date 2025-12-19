import pandas as pd
import subprocess
import sys
import argparse
from concurrent.futures import ProcessPoolExecutor


def fq_sub(raw_path, spname, out_path, sub_num):
    if sub_num < 1000:
        num_reads = round(sub_num*10**9/150/2)
    else
        num_reads = sub_num
    print(num_reads)
    cmd1 = f"mkdir -p {out_path}"
    cmd2 = f"seqtk sample -s 100 {raw_path}/{spname}*R1*gz {num_reads} | gzip > {out_path}/{spname}_sub_{sub_num}G_R1.fq.gz"
    cmd3 = f"seqtk sample -s 100 {raw_path}/{spname}*R2*gz {num_reads} | gzip > {out_path}/{spname}_sub_{sub_num}G_R2.fq.gz"
    try:
        subprocess.check_call(cmd1, shell=True)
        subprocess.check_call(cmd2, shell=True)
        subprocess.check_call(cmd3, shell=True)
    except subprocess.CalledProcessError as e:
        print(f"Command failed with error: {e}")
        sys.exit(1)

#def test(raw_path, spnmae, out_path, sub_num):
#    print("test")

def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, sep='\s+', header=None)
    raw_path_list = df_mapfile.iloc[:,0].tolist()
    spname_list = df_mapfile.iloc[:,1].tolist()
    out_path_list = df_mapfile.iloc[:,2].tolist()
    sub_num_list = df_mapfile.iloc[:,3].tolist()

    return raw_path_list,spname_list,out_path_list,sub_num_list

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mapfile', help='mapfile', required=True)
    # add version
    parser.add_argument('--version', action='version', version='1.0')

    args = parser.parse_args()
    raw_path_list,spname_list,out_path_list,sub_num_list = parse_mapfile(args.mapfile)

    with ProcessPoolExecutor(max_workers = 1) as executor:
        executor.map(fq_sub, raw_path_list,spname_list,out_path_list,sub_num_list)


if __name__ == '__main__':
    main()