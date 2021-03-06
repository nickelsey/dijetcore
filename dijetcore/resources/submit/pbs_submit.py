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

def find_output(configFile):
    jsonfile = open(configFile, 'r')
    json_object = json.load(jsonfile)
    possible_keys = ['outdir', 'output', 'outputDir', 'output_dir', 'output_directory']
    for key in possible_keys:
        if key in json_object:
            return json_object[key]

def generate_submission(executable, input, config, qwrap, execpath, outstream, errstream,
                        name, jobid, prio=0, mem=2, nodes=1, ppn=1, queue='wsuq'):
    submit = 'qsub -V -p {} -lmem={}G -lnodes={}:ppn={} -q {} -o {} -e {} -N {}_{} -- {} {} {} --logtostderr --name={} --id={} --config={} --input={}'.format(
        prio, mem, nodes, ppn, queue, outstream, errstream, name, jobid, qwrap, execpath, executable, name, jobid, config, input)
    return submit


def checkstatus(jobstatus):
    loop = False
    for i in range(len(jobstatus)):
        if jobstatus[i] != 2:
            loop = True
            break
    return loop


def updatestatus(jobstatus, outdir, name):
    from ROOT import TFile
    print("Updating job status")
    print("Total: " + str(len(jobstatus)))

    # get the qstat job listing
    proccommand = 'qstat | grep dx5412'
    proc = subprocess.Popen(proccommand, stdout=subprocess.PIPE, shell=True)
    qstat_result = proc.stdout.read()

    for i in range(len(jobstatus)):

        # if job is completed, we don't need to check again
        if jobstatus[i] == 2:
            continue

        # check if the job is still underway
        jobinprocess = qstat_result.find((name + '_' + str(i)).encode())
        if jobinprocess >= 0:
            jobstatus[i] = 1
            continue

        # if the job is not still underway,
        # check to see if the job has completed properly
        # if not, mark to resubmit
        if outdir.startswith('/'):
            filename = outdir + '/' + name + str(i) + '.root'
        else:
            filename = os.getcwd() + '/' + outdir + '/' + name + str(i) + '.root'

        if os.path.isfile(filename):
            outputfile = TFile(filename, "READ")
            if outputfile.IsZombie():
                print("job " + str(i+1) + " of " + str(len(jobstatus)) +
                      " complete: file is zombie, resubmit")
                jobstatus[i] = 0
                os.remove(filename)
            elif outputfile.IsOpen():
                print("job " + str(i+1) + " of " +
                      str(len(jobstatus)) + " complete: ROOT file healthy")
                print(filename)
                jobstatus[i] = 2
                outputfile.Close()
            else:
                print("job " + str(i+1) + " of " + str(len(jobstatus)) +
                      " undefined file status, resubmit")
                jobstatus[i] = 0
        else:
            print("undefined status: job " + str(i+1) + " of " +
                  str(len(jobstatus)) + " marked for submission")
            jobstatus[i] = 0

    return jobstatus


def main(args):
    # if there are no input files, exit
    files = args.strings
    if files is None:
        return

    # get max number of jobs to be submitted to pbs
    # at a time
    maxjobs = args.maxjobs

    # find binary executable and qwrap script
    execpath = os.getcwd()
    executable = args.binary
    qwrap = args.qwrap

    # find full paths and check they exist
    executable = os.path.abspath(executable)
    qwrap = os.path.abspath(qwrap)

    if not os.path.isfile(executable):
        print('executable: {} does not exist'.format(executable))
        return
    if not os.path.isfile(qwrap):
        print('qwrap: {} does not exist'.format(qwrap))
        return

    # we need to do our own book keeping
    # for when  job is active and when it
    # has completed successfully

    # 0 not submitted, 1 for running, 2 for complete,
    # and -1 for failed - resubmit
    jobstatus = [0 for i in range(len(files))]

    # count the number of qsub submission failures
    qsubfail = 0

    # find our output directory
    output = find_output(args.configFile)

    while checkstatus(jobstatus):

        # update status of all jobs & output
        jobstatus = updatestatus(jobstatus, output, args.name)

        # if we have completed all jobs, exit
        if len(jobstatus) == jobstatus.count(2):
            break

        # find the number of jobs still running via qstat
        # if its at the maximum set jobsactive or jobsqueue,
        # then pause
        jobsactive = jobstatus.count(1)
        while jobsactive >= maxjobs:
            print("reached max number of active jobs: pausing")
            time.sleep(60)
            jobstatus = updatestatus(jobstatus, output, args.name)
            jobsactive = jobstatus.count(1)

        # now submit jobs up to maxjobs - jobsqueued
        njobs = maxjobs - jobsactive

        for i in range(len(jobstatus)):
            if njobs == 0:
                break
            if jobstatus[i] == 1 or jobstatus[i] == 2:
                continue

            # build log & error locations
            outstream = "log/" + args.name + str(i) + ".log"
            errstream = "log/" + args.name + str(i) + ".err"

            qsub = generate_submission(executable, files[i], args.configFile, qwrap, execpath, outstream,
                                       errstream, args.name, i, args.priority, args.mem, args.nodes, 
                                       args.ppn, args.queue)
            print("submitting job: ")
            print(qsub)
            ret = subprocess.Popen(qsub, shell=True)
            ret.wait()
            if ret.returncode == 0:
                jobstatus[i] = 1
            else:
                print("warning: qsub submission failed")
                qsubfail = qsubfail + 1
            njobs = njobs - 1

        if qsubfail > args.max_failures:
            print("qsub failure too many times - exiting")
            return

        # wait one minute before rechecking
        print("finished round of submissions: pausing")
        time.sleep(60)

    print("all jobs completed: exiting")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Submit jobs via PBS & resubmit if necessary')
    parser.add_argument('strings', metavar='S',
                        nargs='*', help=' input files ')
    parser.add_argument('--max_failures', default=200,
                        help=' number of job submissions that can fail before forcing an exit')
    parser.add_argument('--name', default='job_',
                        help=' job name (identifier for qsub)')
    parser.add_argument('--mem', type=int, default=2,
                        help=' memory required per job [GB]')
    parser.add_argument('--nodes', type=int, default=1,
                        help=' number of nodes required per job')
    parser.add_argument('--ppn', type=int, default=1,
                        help=' number of processors per node required per job')
    parser.add_argument('--priority', type=int,
                        default=0, help=' queue priority')
    parser.add_argument('--queue', default='erhiq',
                        help=' queue to submit jobs to')
    parser.add_argument('--maxjobs', type=int, default=100,
                        help=' max number of jobs to have in running or queue states')
    parser.add_argument('--binary', default='', required=True,
                        help='executable with relative path')
    parser.add_argument('--qwrap', default='submit/qwrap.sh',
                        help='qwrap file used during pbs submission')
    parser.add_argument('--configFile', default='config.json',
                        help='JSON configuration file for the specified binary')

    if len(sys.argv) < 2:
        parser.print_usage()
        sys.exit(1)

    main(parser.parse_args())
