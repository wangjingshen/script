import os
import time
import sys

# Define the keywords and time threshold
dir_keywords = [
    ".barcode",
    ".cutadapt",
    ".star",
    ".featureCounts",
    ".assemble",
    ".convert",
    ".summarize",
    ".consensus",
    ".prep_map",
    ".01.Main.anno",
    "seqkit_res",
    "all_matrix",
]
file_suffix = [".bam", ".fq", ".fa",".arrow", "result_mis_2.txt", "annoate_simplify.tsv","release.tar.gz","no_probe_polyT.txt","raw.readnamepath","readname","unmap_cbub.fasta"]
time_threshold_days = 20
time_threshold_seconds = time_threshold_days * 24 * 60 * 60


current_time = time.time()
current_uid = os.getuid()


# Function to check if a file should be deleted
def should_delete_file(file_path):
    try:
        if os.stat(file_path).st_uid != current_uid:
            return False
    except PermissionError:
        return False
    # Check file modification time
    file_mtime = os.path.getmtime(file_path)
    if (current_time - file_mtime) < time_threshold_seconds:
        return False
    # Check if file or directory contains any of the keywords
    for dir_keyword in dir_keywords:
        if dir_keyword in file_path:
            return True
    for suffix in file_suffix:
        if file_path.endswith(suffix):
            return True
    return False


# Function to delete files based on the criteria
def delete_files(base_dirs):
    for base_dir in base_dirs:
        print(f"Scanning {base_dir}")
        for root, dirs, files in os.walk(base_dir):
            for file in files:
                file_path = os.path.join(root, file)
                try:
                    if should_delete_file(file_path):
                        print(f"Deleting: {file_path}")
                        os.remove(file_path)
                except Exception:
                    pass


def main():
    base_directories = sys.argv[1].split(",")
    delete_files(base_directories)
    print("File deletion process completed.")


if __name__ == "__main__":
    main()
