#!/usr/bin/env python

import argparse
import subprocess

def get_promoter_enhancer(gff, outdir):
        cmd1 = (
            f'cat {gff} | grep promoter | '
            r'''awk -F "\t" -v OFS="\t" '{print $1,$4,$5}' | '''
            f'sort -k1,1 -k2,2n > {outdir}/promoter.bed '
        )
        cmd2 = (
            f'cat {gff} | grep enhancer | '
            r'''awk -F "\t" -v OFS="\t" '{print $1,$4,$5}' | '''
            f'sort -k1,1 -k2,2n > {outdir}/enhancer.bed '
        )

        subprocess.check_call(cmd1, shell=True)
        subprocess.check_call(cmd2, shell=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gff', help='gff', required=True)
    parser.add_argument('--outdir', help='outdir', required=True)
    
    args = parser.parse_args()
    get_promoter_enhancer(args.gff, args.outdir)

if __name__ == '__main__':
    main()