#!/usr/bin/env python3

import sys
from os import path, mkdir, getcwd
import subprocess
import shutil
from glob import glob

Reg_cen_list = ['0-5', '5-10', '10-20', '20-30', '30-40', '40-50',
                '50-60', '60-70', '70-80', '80-90', '0-100']
PHOBOS_cen_list = ['0-6', '6-15', '15-25', '25-35', '35-45', '45-55', '45-50']
SPS_cen_list = ['5-12.5', '12.5-23.5', '23.5-33.5', '33.5-43.5']  # SPS PbPb
PHENIX_cen_list = ['20-40', '40-60', '60-88']                     # PHENIX dAu
STAR_cen_list = ['0-10', '10-40', '40-80']                        # STAR v1
centrality_cut_list = (['10-40']
                       #+Reg_cen_list
                       #+ PHOBOS_cen_list
                       #+ SPS_cen_list
                       #+ PHENIX_cen_list
                       #+ STAR_cen_list
                      )

def generate_script(folder_name, ecm, centrality_label, sys_label, nev):
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
##PBS -q ecsxq
#PBS -q wsuq
##PBS -m bea
##PBS -M chunshen1987@gmail.com

cd {2:s}
rm -fr data
mkdir data
python3 generateEbeprofiles.py -ecm {3:g} -cen {4:s} -collision_system {5:s} -nev {6:d}

""".format(working_folder.split('/')[-1], walltime, working_folder, ecm,
           centrality_label, sys_label, nev))
    script.close()


def generate_event_folder(working_folder, event_id, ecm, collision_system, nev):
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
    shutil.copyfile(
        path.join(code_path, 'scripts/generateEbeprofiles.py'),
        path.join(working_folder, 'generateEbeprofiles.py'))
    subprocess.call("ln -s {0:s} {1:s}".format(
        path.join(code_path, "scripts/centrality_cut_tables"),
        path.join(working_folder, "centrality_cut_tables")), shell=True)
    generate_script(working_folder, ecm, centrality_cut_list[event_id],
                    collision_system, nev)

if __name__ == "__main__":
    try:
        folder_name = str(sys.argv[1])
        ecm = float(sys.argv[2])
        collision_system = str(sys.argv[3])
        nev = int(sys.argv[4])
    except IndexError:
        print("{} working_folder ecm collision_system nev".format(
                                                        str(sys.argv[0])))
        exit(0)

    ncore = len(centrality_cut_list)
    mkdir(folder_name)
    for icore in range(ncore):
        generate_event_folder(folder_name, icore, ecm, collision_system, nev)

