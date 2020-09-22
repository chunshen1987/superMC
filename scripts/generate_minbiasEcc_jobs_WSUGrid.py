#!/usr/bin/env python3

import sys
from os import path, mkdir, getcwd
import subprocess
import shutil
from glob import glob

def generate_script(folder_name, ecm):
    working_folder = path.join(path.abspath('./'), folder_name)
    walltime = '20:00:00'

    script = open(path.join(working_folder, "submit_job.pbs"), "w")
    script.write(
"""#!/usr/bin/env bash
#SBATCH --job-name {0:s}
#SBATCH -q primary
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=4G
#SBATCH --constraint=intel
#SBATCH -t {1:s}
#SBATCH -e job.err
#SBATCH -o job.log

cd {2:s}
rm -fr data
mkdir data
./superMC.e ecm={3:g} operation=9 ecc_from_order=1 ecc_to_order=9 bmin=0. bmax=20. finalFactor=1.0

""".format(working_folder.split('/')[-1], walltime, working_folder, ecm))
    script.close()


def generate_event_folder(working_folder, event_id, ecm):
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
    generate_script(working_folder, ecm)

if __name__ == "__main__":
    try:
        folder_name = str(sys.argv[1])
        ncore = int(sys.argv[2])
        ecm = float(sys.argv[3])
    except IndexError:
        print("{} working_folder num_of_cores ecm".format(str(sys.argv[0])))
        exit(0)

    mkdir(folder_name)
    for icore in range(ncore):
        generate_event_folder(folder_name, icore, ecm)

