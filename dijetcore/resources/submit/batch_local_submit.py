# pbs_submit.py

# convenience script to batch submit and monitor
# job status, automatically resubmitting jobs that
# do not produce the expected output files

from __future__ import print_function, division

import sys
import os
import argparse
import subprocess
import time
import json
import multiprocessing
import atexit

def find_output(configFile):
    jsonfile = open(configFile, 'r')
    json_object = json.load(jsonfile)
    possible_keys = ['outdir', 'output', 'outputDir', 'output_dir', 'output_directory']
    for key in possible_keys:
        if key in json_object:
            return json_object[key]
    return None

def generate_submission(executable, name, index, config):
    submit = [executable, '-name='+name, '-id='+str(index), '-config='+config]
    return submit

def check_status(job_status):
    return job_status.count(1) != len(job_status)

def update_status(job_status, threads):
    for i in range(len(threads)):
        t = threads[i]
        if t is not None:
            job_status[i] = 0
            ret_val = t.poll()
            if ret_val is not None:
                if ret_val is not 0:
                    print('job {} exited with non-zero status - possible error'.format(i))
                threads[i] = None
                job_status[i] = 1

def cleanup(threads):
    for t in threads:
        if t is not None:
            t.kill()

def main(args):
    # if there are no input files, exit
    files = args.strings
 
    # get number of jobs - its either defined by the number of files
    # or if there are no input files, by the numJobs arg
    num_jobs = 0
    if files:
        num_jobs = len(files)
    else:
        num_jobs = args.numJobs

    # get max number of jobs to be submitted at once
    # if the number is larger than the number of cores, restrict
    max_jobs = args.maxConcurrentJobs
    if max_jobs > multiprocessing.cpu_count():
        max_jobs = multiprocessing.cpu_count();

    # find binary executable 
    exec_path = os.getcwd()
    executable = args.binary

    # find full paths and check they exist
    executable = os.path.abspath(executable)

    if not os.path.isfile(executable):
        print('executable: {} does not exist'.format(executable))
        return

    # find our output directory
    output = find_output(args.configFile)

    if output:
        try:
            os.makedirs(output, exist_ok=True)
        except OSError:
            print('could not create output directory: {}, exiting'.format(output))
            return

    # find our log directory
    log_dir = args.logDir

    if log_dir:
        try:
            os.makedirs(log_dir, exist_ok=True)
        except OSError:
            print('could not create log directory: {}, exiting'.format(log_dir))
            return

    job_status = [-1 for i in range(num_jobs)]
    threads = [None for i in range(num_jobs)]
    current_job_id = 1

    # in case of unannounced exit - cleanup  
    atexit.register(cleanup, threads)

    # start while loop, continue until all jobs have successfully completed
    while check_status(job_status):
        
        update_status(job_status, threads)
        
        jobs_active = job_status.count(0)
        
        while jobs_active >= max_jobs:
            print("reached max number of active jobs: pausing")
            time.sleep(30)
            update_status(job_status, threads)
            jobs_active = job_status.count(0)

        # now submit jobs up to max_jobs - jobs_active or the number of jobs left
        # whichever is smaller
        n_jobs_submit = max_jobs - jobs_active
        
        # start submissions
        for i in range(num_jobs):
            if n_jobs_submit <= 0:
                break
            if job_status[i] > -1:
                continue
            
            # this index will be submitted - create log directory 
            std_out = open('{}/{}{}.log'.format(log_dir, args.name, current_job_id + args.startIndex), 'w')
            std_err = open('{}/{}{}.err'.format(log_dir, args.name, current_job_id + args.startIndex), 'w')

            # create job submission for Popen
            job_submit = generate_submission(executable, args.name, current_job_id + args.startIndex, args.configFile)
            
            print('submitting job: {}\nlog directory: {}\nerr directory: {}'.format(job_submit, std_out.name, std_err.name))
            # create subprocess and store
            threads[i] = subprocess.Popen(job_submit, stdout=std_out, stderr=std_err)

            # increase current job id
            current_job_id += 1
            n_jobs_submit -= 1

        update_status(job_status, threads)
        
        # wait 30 seconds before rechecking
        print("finished round of submissions: pausing")
        time.sleep(30)

    print('all jobs completed: exiting')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Run batches of jobs locally')
    parser.add_argument('strings', metavar='S',
                        nargs='*', help=' input files ')
    parser.add_argument('--maxConcurrentJobs', type=int, default=10,
                        help=' number of job submissions to run concurrently')
    parser.add_argument('--name', default='job_',
                        help=' job name (identifier for output)')
    parser.add_argument('--startIndex', type=int, default=0,
                        help='index to start submission at')
    parser.add_argument('--numJobs', type=int, default=1,
                        help=' number of times to submit the job: not used if input files are specified')
    parser.add_argument('--binary', default='', required=True,
                        help='executable with relative path')
    parser.add_argument('--logDir', default='log', help= 'directory to store logs in')
    parser.add_argument('--configFile', default='config.json',
                        help='JSON configuration file for the specified binary')

    if len(sys.argv) < 2:
        parser.print_usage()
        sys.exit(1)

    main(parser.parse_args())
