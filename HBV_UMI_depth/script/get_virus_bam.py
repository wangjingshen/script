import sys
import pysam
import pandas as pd
from collections import defaultdict
import argparse
import os
import glob
import subprocess
from concurrent.futures import ProcessPoolExecutor

def get_virus_barcode_umi(virus_tsne, virus_read_count, otsu_min_support_read):
    df_virus_tsne = pd.read_csv(virus_tsne, sep="\t", index_col=0).fillna(0)
    virus_barcode = list(df_virus_tsne['barcode'][df_virus_tsne.UMI!=0])
    df_virus_read_count = pd.read_csv(virus_read_count, sep="\t")
    df_virus_read_count = df_virus_read_count[(df_virus_read_count.read_count >= otsu_min_support_read) & (df_virus_read_count.barcode.isin(virus_barcode)) ]
    virus_barcode_umi = list(df_virus_read_count.barcode + "_" + df_virus_read_count.UMI)
    return(virus_barcode_umi)


class Bam_virus:
    def __init__(self, celescope_path, otsu_min_support_read, outdir, name):
        self.input_bam = glob.glob(f"{celescope_path}/*star*/*Aligned.sortedByCoord.out.bam")[0]
        self.virus_tsne = glob.glob(f"{celescope_path}/*analysis_capture_virus*/*virus_tsne.tsv")[0]
        self.virus_read_count = glob.glob(f"{celescope_path}/*count_capture_virus/*virus_read_count.tsv")[0]
        self.otsu_min_support_read = int(otsu_min_support_read)
        self.outdir = outdir
        self.name = name
        if not os.path.exists(self.outdir):
            os.system(f"mkdir -p {self.outdir}")

    def virus_bam(self):
        virus_barcode_umi = get_virus_barcode_umi(self.virus_tsne, self.virus_read_count, self.otsu_min_support_read)
        input_bam = pysam.AlignmentFile(self.input_bam, "rb")
        out_bam = pysam.AlignmentFile(f'{self.outdir}/{self.name}_virus.bam', "wb", template = input_bam)
        ref_record_dict = defaultdict(list)
        for i in input_bam:
            ref_record_dict["qname"].append(i.qname)
            ref_record_dict["record"].append(i)
        ref_record_df = pd.DataFrame(ref_record_dict)
        ref_record_df['barcode_umi'] = ref_record_df['qname'].str.split("_").str[0] + "_" + ref_record_df['qname'].str.split("_").str[1]
        ref_record_df_filter = ref_record_df[ref_record_df.barcode_umi.isin(virus_barcode_umi)]
        for _, row in ref_record_df_filter.iterrows():
            out_bam.write(row["record"])
        input_bam.close()
        out_bam.close()
            
    def bam2tsv(self):
        cmd = (f'samtools view {self.outdir}/{self.name}_virus.bam | '
               r'''awk -F "\t" '{print $1}' | awk -F "_" '{print $1"_"$2"\t"$3}' > '''
               f'{self.outdir}/{self.name}_virus_id.tsv'
              )
        subprocess.check_call(cmd, shell = True)
    
    def run(self):
        self.virus_bam()
        self.bam2tsv()

def parse_mapfile(mapfile):
    df_mapfile = pd.read_csv(mapfile, sep='\t')
    celescope_path_list = df_mapfile['celescope_path']
    otsu_min_support_read_list = df_mapfile['otsu_min_support_read']
    outdir_list = df_mapfile['outdir']
    name_list = df_mapfile['name']
    
    return celescope_path_list, otsu_min_support_read_list, outdir_list, name_list

def run_single(celescope_path, otsu_min_support_read, outdir, name):
    runner = Bam_virus(celescope_path, otsu_min_support_read, outdir, name)
    runner.run()

def main():
    mapfile = sys.argv[1]
    celescope_path_list, otsu_min_support_read_list, outdir_list, name_list = parse_mapfile(mapfile)
    with ProcessPoolExecutor(max_workers = 1) as executor:
        executor.map(run_single, celescope_path_list, otsu_min_support_read_list, outdir_list, name_list)


if __name__ == '__main__':
    main()