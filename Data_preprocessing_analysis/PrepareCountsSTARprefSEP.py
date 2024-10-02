
import os
import subprocess
import time
import sys
import pandas as pd
import shutil
import psutil

"""
Example commands:

conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis
/home/aleksandr/miniconda3/envs/tf-py38/bin/python PrepareCountsSTARprefSEP.py test_id.txt PAIRED GSE102556_


"""

def kill_process_and_children(pid):
    try:
        parent = psutil.Process(pid)
        children = parent.children(recursive=True)
        for child in children:
            child.kill()
        parent.kill()
    except psutil.NoSuchProcess:
        pass

# /home/aleksandr/miniconda3
# $CONDA_PREFIX/etc/profile.d/conda.sh && conda activate test4 && conda env list

#### Defining function to get file ####
def try_to_prefetch(ID, folder, timeout_par=1200):

    abs_path = os.path.abspath("")

    Error_counter_connection = 0
    Error_flag = False
    TRY_reprefetch_flag = False

    request_template = f"""
    source /home/aleksandr/miniconda3/etc/profile.d/conda.sh && conda activate tf-py38 ;
    source /etc/profile.d/sra-tools.sh ;
    cd {abs_path} ;
    prefetch {ID} --output-directory {folder} --max-size 60g
    """

    print(request_template)

    try: 

        os.makedirs(folder, exist_ok=True)
        process = subprocess.Popen(request_template, shell=True, executable="/bin/bash")
        process.wait(timeout=timeout_par)
        print(f"Download was successful")

    except subprocess.CalledProcessError as error:
        print(error)
        kill_process_and_children(process.pid)
        Error_counter_connection += 1
        Error_flag = True
        TRY_reprefetch_flag = True

    except Exception as e:
        print(f"An error occurred: {str(e)}")
        kill_process_and_children(process.pid)
        shutil.rmtree(folder)
        Error_counter_connection += 1
        Error_flag = True
        TRY_reprefetch_flag = True

    finally:
        kill_process_and_children(process.pid)
        if process.stdout:
            process.stdout.close()
        if process.stderr:
            process.stderr.close()
        if process.stdin:
            process.stdin.close()

        del process

    while Error_flag:

        print(Error_counter_connection," FAILED_TO_DOWNLOAD_FROM ", ID, " ", time.strftime("%H:%M:%S",time.localtime()))
        time.sleep(30 + 2**Error_counter_connection)

        try:

            if TRY_reprefetch_flag:

                print(f"Attempting to prefetch again without disruption")
                process = subprocess.Popen(request_template, shell=True, executable="/bin/bash")
                process.wait(timeout=timeout_par)
                print(f"Download was successful")
                Error_flag = False
                TRY_reprefetch_flag = False

            else:
                os.makedirs(folder, exist_ok=True)
                process = subprocess.Popen(request_template, shell=True, executable="/bin/bash")
                process.wait(timeout=timeout_par)
                print(f"Download was successful")
                Error_flag = False

        except subprocess.CalledProcessError as error:
            print(error)
            kill_process_and_children(process.pid)
            shutil.rmtree(folder)
            Error_flag = True
            TRY_reprefetch_flag = False # reprefetch happens only the 1st time!
            Error_counter_connection += 1

        except Exception as e:
            print(f"An error occurred: {str(e)}")
            kill_process_and_children(process.pid)
            shutil.rmtree(folder)
            Error_counter_connection += 1
            Error_flag = True
            TRY_reprefetch_flag = False # reprefetch happens only the 1st time!

        finally:
            kill_process_and_children(process.pid)
            if process.stdout:
                process.stdout.close()
            if process.stderr:
                process.stderr.close()
            if process.stdin:
                process.stdin.close()
            del process



def try_to_fasterq_dump(ID, folder, timeout_par=2000):

    abs_path = os.path.abspath("")

    Error_counter_connection = 0
    Error_flag = False

    request_template = f"""
    source /home/aleksandr/miniconda3/etc/profile.d/conda.sh && conda activate tf-py38 ;
    source /etc/profile.d/sra-tools.sh ;
    cd {abs_path} ;
    cd {folder} ;
    fasterq-dump {ID} --threads 9 --force
    """

    print(request_template)

    try: 

        os.makedirs(folder, exist_ok=True)
        process = subprocess.Popen(request_template, shell=True, executable="/bin/bash")
        process.wait(timeout=timeout_par)
        print(f"FASTQ generated successfully")

    except subprocess.CalledProcessError as error:
        print(error)
        kill_process_and_children(process.pid)
        Error_counter_connection += 1
        Error_flag = True

    except Exception as e:
        print(f"An error occurred: {str(e)}")
        kill_process_and_children(process.pid)
        Error_counter_connection += 1
        Error_flag = True

    finally:
        kill_process_and_children(process.pid)
        if process.stdout:
            process.stdout.close()
        if process.stderr:
            process.stderr.close()
        if process.stdin:
            process.stdin.close()

        del process

    while Error_flag:

        print(Error_counter_connection," FAILED_TO_GENERATE_FASTQ:  ", ID, " ", time.strftime("%H:%M:%S",time.localtime()))
        time.sleep(30 + 2**Error_counter_connection)

        try:
            os.makedirs(folder, exist_ok=True)
            process = subprocess.Popen(request_template, shell=True, executable="/bin/bash")
            process.wait(timeout=timeout_par)
            print(f"FASTQ generated successfully")
            Error_flag = False

        except subprocess.CalledProcessError as error:
            print(error)
            kill_process_and_children(process.pid)
            Error_flag = True
            Error_counter_connection += 1

        except Exception as e:
            print(f"An error occurred: {str(e)}")
            kill_process_and_children(process.pid)
            Error_counter_connection += 1
            Error_flag = True

        finally:
            kill_process_and_children(process.pid)
            if process.stdout:
                process.stdout.close()
            if process.stderr:
                process.stderr.close()
            if process.stdin:
                process.stdin.close()
            del process


def read_file(file_path):

    file_lines = []

    with open(file_path, "r") as file:
        for line in file:
            file_lines.append(line.strip())

    return file_lines

def save_file(file_path, object):
    with open(file_path, "w") as file:
        for line in object:
            file.write(line + "\n")

def list_files(directory):
    paths = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            paths.append(os.path.join(root, file))
    return paths

def run_STAR_request(fastq_folder, STAR_dir, paired):

    abs_path = os.path.abspath("")

    files_fastq = list_files(fastq_folder)
    files_fastq = [x for x in files_fastq if ".fastq" in x]
    files_fastq.sort()

    if paired:

        STAR_request = f"""
        source /home/aleksandr/miniconda3/etc/profile.d/conda.sh && conda activate tf-py38 ;
        cd {abs_path} ;
        STAR --runMode alignReads --runThreadN 9 --genomeDir genome_index_STAR_hg19/ --outSAMtype BAM SortedByCoordinate \
        --readFilesIn {files_fastq[0]} {files_fastq[1]} \
        --outFileNamePrefix {STAR_dir}/ --quantMode GeneCounts --sjdbGTFfile genome_files_hg19/Homo_sapiens.GRCh37.87.gtf ;

        """

    else:

        STAR_request = f"""
        source /home/aleksandr/miniconda3/etc/profile.d/conda.sh && conda activate tf-py38 ;
        cd {abs_path} ;
        STAR --runMode alignReads --runThreadN 9 --genomeDir genome_index_STAR_hg19/ --outSAMtype BAM SortedByCoordinate \
        --readFilesIn {files_fastq[0]} \
        --outFileNamePrefix {STAR_dir}/ --quantMode GeneCounts --sjdbGTFfile genome_files_hg19/Homo_sapiens.GRCh37.87.gtf ;

        """

    print(STAR_request)
    process = subprocess.Popen(STAR_request, shell=True, executable="/bin/bash")
    process.wait()
    kill_process_and_children(process.pid)
    if process.stdout:
        process.stdout.close()
    if process.stderr:
        process.stderr.close()
    if process.stdin:
        process.stdin.close()
    del process

#### Parsing arguments ####
ARG_LIST = sys.argv

# Run parameters
# 1st ID txt.file
# 2nd Paired PAIRED or SINGLE
# 3rd Directory prefix (For example: "GSE45345_")

ID_ARG = str(ARG_LIST[1])

PAIRED_ARG = str(ARG_LIST[2])
if PAIRED_ARG == "PAIRED":
    PAIRED_ARG = True
else:
    PAIRED_ARG = False

DIR_ARG = str(ARG_LIST[3])


#### Downloading reads with fasterq-dump ####

ids = read_file(ID_ARG)

fastq_dir = DIR_ARG + "fastq"

STAR_dir = DIR_ARG + "STAR"

OUTPUT_dir = DIR_ARG + "OUTPUT"
os.makedirs(OUTPUT_dir, exist_ok=True)

for ID in ids:

    print("---------------------")
    print(f"Working on: {ID}")
    print(time.strftime("%H:%M:%S",time.localtime()))

    os.makedirs(STAR_dir, exist_ok=True)

    # run preprocessing things
    try_to_prefetch(ID, fastq_dir)
    try_to_fasterq_dump(ID, fastq_dir)
    run_STAR_request(fastq_dir, STAR_dir, PAIRED_ARG)

    # move files to desired location
    STAR_files = list_files(STAR_dir)
    STAR_counts = [x for x in STAR_files if "ReadsPerGene.out.tab" in x]
    STAR_logs = [x for x in STAR_files if "Log.final.out" in x]

    logs_list = read_file(STAR_logs[0])
    counts_df = pd.read_csv(STAR_counts[0], sep="\t", header=None)
    counts_df.columns = ["gene_ID", "Unstranded", "Strand_1", "Strand_2"]

    log_path = os.path.join(OUTPUT_dir, ID + "_star_logs.txt")
    counts_path = os.path.join(OUTPUT_dir, ID + "_star_counts.csv")

    save_file(log_path, logs_list)
    counts_df.to_csv(counts_path)

    shutil.rmtree(fastq_dir)
    shutil.rmtree(STAR_dir)

    print(f"Working on: {ID} is DONE!")
    print(time.strftime("%H:%M:%S",time.localtime()))

