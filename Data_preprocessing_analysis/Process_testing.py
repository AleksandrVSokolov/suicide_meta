import subprocess
import select
import sys
import time

"""
Example commands:

conda activate tf-py38
cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis
/home/aleksandr/miniconda3/envs/tf-py38/bin/python Process_testing.py


"""


def stop_process_on_output(command, target_output, interval=10):
    # Start the process
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, shell=True, executable="/bin/bash")
    
    try:

        timer = 0
        timer_stop = 1200

        while True:
            # Use select to wait for the process output with a timeout

            line = process.stdout.readline()
            if line == '' and process.poll() is not None:
                print("Process completed before timeout.")
                break  # Process has completed, so exit the loo
            
            print(line, end='')  # Print the output (optional)

            # Check if the target output is in the current line
            if any(target_output in line for target_output in patterns):

                print("Target output found. Terminating the process.")
                process.terminate()  # Terminate the process
                break

            timer += 1
            time.sleep(1)

            if timer > timer_stop:
                print("Timeout is reached -> terminating the process")
                process.terminate()  # Terminate the process
                break


    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)

    finally:
        # Ensure the process is terminated and cleaned up
        process.stdout.close()
        process.wait()

command = """
    source /home/aleksandr/miniconda3/etc/profile.d/conda.sh && conda activate tf-py38 ;
    source /etc/profile.d/sra-tools.sh ;
    cd /home/aleksandr/Desktop/WORK/OLINK_suicide_PSY_project/Data_preprocessing_analysis ;
    mkdir GSE247998_pref ;
    prefetch SRR26854272 --output-directory GSE247998_pref --max-size 60g ;
    cd GSE247998_pref ;
    fasterq-dump SRR26854272 --threads 9 --force
    """  
patterns = ["HTTPS", "Error", "ERROR", "fail"]  # The text to match
stop_process_on_output(command, patterns)