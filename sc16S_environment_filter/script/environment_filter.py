import pandas as pd
import argparse
import subprocess


def environment_filter(environment_genus, sample):
    raw_matrix_file = f"{sample}/outs/{sample}_raw_UMI_matrix.tsv.gz"
    filter_matrix_file = f"{sample}/outs/{sample}_filter_UMI_matrix.tsv.gz"

    df = pd.read_csv(raw_matrix_file, sep="\t", index_col=0)
    environment_genus = pd.read_csv(environment_genus, header = None)
    environment_genus = environment_genus.iloc[:,0]

    df_filter = df[~df.index.isin(environment_genus)]  # rm environment
    df_filter.to_csv(filter_matrix_file, sep="\t")


if __name__ == '__main__':
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--environment_genus', help="environment_genus", required=True)
    parser.add_argument("--sample", help="sample", required=True)

    # add version
    parser.add_argument('--version', action='version', version='1.0')
    args = parser.parse_args()

    environment_filter(args.environment_genus, args.sample)