#!/usr/bin/env python3

import sys
from os import path, mkdir, getcwd
import subprocess
import shutil
from glob import glob

def generate_script(folder_name):
    working_folder = path.join(path.abspath('./'), folder_name)
    walltime = '10:00:00'

    script = open(path.join(working_folder, "submit_job.pbs"), "w")
    script.write(
"""#!/usr/bin/env bash
#PBS -N {0:s}
#PBS -l select=1:ncpus=1:mem=4GB:cpu_type=Intel
##PBS -l walltime={1:s}
#PBS -S /bin/bash
#PBS -e test.err
#PBS -o test.log
#PBS -q ecsxq
##PBS -m bea
##PBS -M chunshen1987@gmail.com

cd {2:s}
rm -fr data
mkdir data
./superMC.e

""".format(working_folder.split('/')[-1], walltime, working_folder))
    script.close()


def generate_event_folder(working_folder, event_id):
    working_folder = path.join(
            path.abspath(working_folder), "superMC_{}".format(event_id))
    mkdir(working_folder)
    code_path = path.join(getcwd(), "codes/superMC")
    subprocess.call("ln -s {0:s} {1:s}".format(
        path.join(code_path, "superMC.e"),
        path.join(working_folder, "superMC.e")), shell=True)
    subprocess.call("ln -s {0:s} {1:s}".format(
        path.join(code_path, "EOS"),
        path.join(working_folder, "EOS")), shell=True)
    subprocess.call("ln -s {0:s} {1:s}".format(
        path.join(code_path, "tables"),
        path.join(working_folder, "tables")), shell=True)
    shutil.copyfile(
        path.join(code_path, 'parameters.dat'),
        path.join(working_folder, 'parameters.dat'))
    generate_script(working_folder)

if __name__ == "__main__":
    try:
        folder_name = str(sys.argv[1])
        ncore = int(sys.argv[2])
    except IOError:
        print("./generate_jobs_guillimin.py working_folder num_of_cores")
        exit(0)

    for icore in range(ncore):
        generate_event_folder(folder_name, icore)

