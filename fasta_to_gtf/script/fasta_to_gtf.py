#!/usr/bin/env python

import os
import sys
import argparse
import subprocess
import logging
from pathlib import Path
from contextlib import contextmanager
import time

dev_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(dev_root))
from utils.utils import mkdir, logger, execute_cmd

# Context manager: temporarily switch the working directory
@contextmanager
def cwd(path: Path):
    origin = Path.cwd()
    try:
        os.chdir(path)
        yield
    finally:
        os.chdir(origin)

class FastaToGtf:
    def __init__(self, genome:Path, fasta:Path, prefix:str, star_path:Path, threads:int=4, force:bool=False):
        self.genome = genome.resolve()
        self.fasta = fasta.resolve()
        self.prefix = prefix
        self.star_path = star_path
        self.threads = threads
        self.force = force
        self.star_out = Path("star")
        self.gtf_out = Path(f'{prefix}.gtf')

    def star_align(self) -> None:
        """STAR to get BAM"""
        bam_file = self.star_out / f"{self.prefix}Aligned.sortedByCoord.out.bam"
        if bam_file.exists() and not self.force:
            logger.info(f"BAM already exists: {bam_file}  (--force to overwrite)")
            return

        mkdir(self.star_out)
        cmd = (
            f"{self.star_path} "
            f"--runThreadN {self.threads} "
            f"--genomeDir {self.genome} "
            f"--readFilesIn {self.fasta} "
            f"--outSAMtype BAM SortedByCoordinate "
            f"--outFileNamePrefix {self.prefix}"
        )
        logger.info("Running STAR...")
        with cwd(self.star_out):
            execute_cmd(cmd)

    def bam_to_bed(self) -> Path:
        """BAM → bed12"""
        bed = Path(f"{self.prefix}.bed")
        bam = f'{self.star_out}/{self.prefix}Aligned.sortedByCoord.out.bam'
        cmd = f"bedtools bamtobed -i {bam} -bed12 > {bed}"
        logger.info("BAM → BED...")
        execute_cmd(cmd)
        return bed

    def bed_to_gtf(self, bed: Path) -> Path:
        """bed → genePred → GTF"""
        gtf = Path(f"{self.prefix}.gtf")
        cmd = (
            f"bedToGenePred {bed} /dev/stdout | "
            f"genePredToGtf file /dev/stdin tmp && "
            f"sed -i 's|/dev/stdin|{self.prefix}_genome|g' tmp && "
            f"awk -F'\\t' -v OFS='\\t' '{{if ($3==\"transcript\") $3=\"gene\"; print}}' tmp > {gtf} && "
            f"rm -f tmp"
        )
        logger.info("BED → GTF...")
        execute_cmd(cmd)
        return gtf

    def run(self) -> None:
        t0 = time.time()
        self.star_align()
        bed = self.bam_to_bed()
        gtf = self.bed_to_gtf(bed)
        logger.info(f"Done → {gtf}  ({time.time() - t0:.1f}s)")


def main():
    parser = argparse.ArgumentParser(description="FASTA → STAR → BAM → BED → GTF")
    sub = parser.add_subparsers(dest="command", required=True)

    p_run = sub.add_parser("run", help="run pipeline")
    p_run.add_argument("--genome", required=True, type=Path, help="genome path")
    p_run.add_argument("--fasta",  required=True, type=Path, help="fasta")
    p_run.add_argument("--prefix", required=True, help="prefix")
    p_run.add_argument("--star",   default="STAR", help="STAR path, /SGRNJ/Public/Software/conda_env/celescope2.1.0/bin/STAR-avx2,\
                                                                    /Public/Software/miniconda2/envs/old/bin/STAR")
    p_run.add_argument("--threads", default=4, type=int, help="threads")
    p_run.add_argument("--force",  action="store_true", help="force get bam")
    args = parser.parse_args()

    if args.command == "run":
        star_path = Path(args.star).expanduser()
        if not star_path.is_file():
            raise FileNotFoundError(f"STAR not found: {star_path}")
        runner = FastaToGtf(args.genome, args.fasta, args.prefix, star_path, args.threads, args.force)
        runner.run()

if __name__ == "__main__":
    main()