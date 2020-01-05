#! /usr/bin/env python

from sys import argv, exit
from os import path
import h5py
import numpy as np

def print_usage():
    """This function prints out help messages"""
    print("Usage: {} ".format(argv[0]) + "database.h5")


def output_centrality_cut_table_iebe_package(database, cut_type,
                                             folder_path, database_name,
                                             alpha=0.0):
    """
    This function output centrality cut table for iEBE package
    :param cut_type: centrality cut type
    """
    if cut_type == "GlauberMixed":
        print(
            'cutting centralities according to {} with alpha = {}....'.format(
                                                            cut_type, alpha)
        )
        centrality_output = (
            open('{0}/iebe_centralityCut_{1}_{2}_alpha_{3}.dat'.format(
                            folder_path, cut_type, database_name, alpha), 'w')
        )
    else:
        print('cutting centralities according to {} ....'.format(cut_type))
        centrality_output = (
            open('{0}/iebe_centralityCut_{1}_{2}.dat'.format(
                            folder_path, cut_type, database_name), 'w')
        )

    coll_data = np.array(hf.get("collision_data"))
    nevent, ncol = coll_data.shape

    sorted_idx = range(nevent)
    if cut_type == 'b':
        sorted_idx = np.argsort(coll_data[:, 0])
        centrality_output.write("#centrality b(fm)\n")
        centrality_output.write("%6.4e %18.8e\n" % (0.0, 0.0))
    elif cut_type == "Npart":
        sorted_idx = np.argsort(-coll_data[:, 1])
        centrality_output.write("#centrality Npart b_min(fm) b_max(fm)\n")
        centrality_output.write("%6.4e  %18.8e  %18.8e  %18.8e\n"
            % (0.0, 500, 0.0, 0.0))
    elif cut_type == "total_entropy":
        sorted_idx = np.argsort(-coll_data[:, 3])
        centrality_output.write(
            "#centrality dS/dy Npart_min Npart_max b_min(fm) b_max(fm)\n")
        centrality_output.write(
            "%6.4e  %18.8e  %18.8e  %18.8e  %18.8e  %18.8e\n"
            % (0.0, 1000000, 500, 500, 0.0, 0.0))
    elif cut_type == "GlauberMixed":
        sorted_idx = np.argsort(
            -((1. - alpha)/2.*coll_data[:, 1] + alpha*coll_data[:, 2])
        )
        centrality_output.write(
            "#centrality dS/dy Npart_min Npart_max b_min(fm) b_max(fm)\n")
        centrality_output.write(
            "%6.4e  %18.8e  %18.8e  %18.8e  %18.8e  %18.8e\n"
            % (0.0, 1000000, 500, 500, 0.0, 0.0))
    else:
        print("invalid cutType", cut_type)
        exit(0)

    centrality_bound = (
        [x*0.1 for x in range(0, 10)] + [y for y in range(1, 101)]
    )
    for icen in range(1, len(centrality_bound)):
        lowerbound = centrality_bound[icen-1]
        upperbound = centrality_bound[icen]
        nsample = int(nevent * (upperbound - lowerbound) / 100) - 1
        noffset = int(nevent * lowerbound / 100)

        fetched_data = coll_data[sorted_idx[noffset:noffset+nsample], :]

        npart_min = min(fetched_data[:, 1])
        npart_max = max(fetched_data[:, 1])
        b_min = min(fetched_data[:, 0])
        b_max = max(fetched_data[:, 0])
        dsdy_min = min(fetched_data[:, 3])
        if cut_type == 'total_entropy':
            centrality_output.write(
                "%6.4e  %18.8e  %18.8e  %18.8e  %18.8e  %18.8e\n"
                % (upperbound, dsdy_min, npart_min, npart_max,
                   b_min, b_max)
            )
        elif cut_type == 'GlauberMixed':
            dsdy_min = min((1. - alpha)/2.*fetched_data[:, 1]
                           + alpha*fetched_data[:, 2])
            centrality_output.write(
                "%6.4e  %18.8e  %18.8e  %18.8e  %18.8e  %18.8e\n"
                % (upperbound, dsdy_min, npart_min, npart_max,
                   b_min, b_max)
            )
        elif cut_type == 'Npart':
            centrality_output.write("%6.4e  %18.8e  %18.8e  %18.8e\n"
                % (upperbound, npart_min, b_min, b_max)
            )
        elif cut_type == 'b':
            centrality_output.write("%6.4e  %18.8e\n"
                % (upperbound, b_max)
            )
    centrality_output.close()


try:
    database_path = path.abspath(str(argv[1]))
    database_name = database_path.split("/")[-1].split(".h5")[0]
    database_directory = "/".join(database_path.split("/")[0:-1])
except IndexError:
    print_usage()
    exit(1)

print("generating centrality cut table from {} ...".format(database_path))
hf = h5py.File(database_path, "r")
output_centrality_cut_table_iebe_package(hf, "total_entropy",
                                         database_directory, database_name)
output_centrality_cut_table_iebe_package(hf, "Npart",
                                         database_directory, database_name)
#output_centrality_cut_table_iebe_package(hf, "GlauberMixed",
#                                         database_directory, database_name,
#                                         0.12)
hf.close()
print("finished.")
