#!/usr/bin/env python3

import numpy as np
import h5py
from os import path
import sys

def print_usage():
    """This function prints out help messages"""
    print("Usage: {} ".format(sys.argv[0]) + "result_folder")

try:
    results_folder = path.abspath(str(sys.argv[1]))
    results_name = results_folder.split("/")[-1]
except IndexError:
    print_usage()
    exit(1)

print("collecting data from {} ...".format(results_folder))
hf = h5py.File("{0}/{1}.h5".format(results_folder, results_name), "w")

print("loading ascii data from disk ...")
data_sn = np.loadtxt("{0}/sn_ecc_eccp_10.dat".format(results_folder))
data_en = np.loadtxt("{0}/en_ecc_eccp_10.dat".format(results_folder))

nev, ncol = data_sn.shape
print("read in {} events.".format(nev))
coll_data = np.zeros([nev, 4])
coll_data[:, 0] = data_sn[:, 39]  # b [fm]
coll_data[:, 1] = data_sn[:, 36]  # Npart
coll_data[:, 2] = data_sn[:, 37]  # Ncoll
coll_data[:, 3] = data_sn[:, 38]  # dS/dy

# save collision data
header_text = "# b[fm]  Npart  Ncoll  dS/dy"
h5data = hf.create_dataset("collision_data", data=coll_data.astype('float32'),
                           compression="gzip", compression_opts=9)
h5data.attrs.create("header", np.string_(header_text))

# save eccentricity data
ecc_data = np.zeros([nev, 18])
for iorder in range(1, 10):
    ecc_data[:, 2*iorder-2] = data_sn[:, 4*iorder-2]
    ecc_data[:, 2*iorder-1] = data_sn[:, 4*iorder-1]
header_text = "# ecc_n_real  ecc_n_imag (n = 1-9)"
h5data = hf.create_dataset("ecc_data_sn", data=ecc_data.astype('float32'),
                           compression="gzip", compression_opts=9)
h5data.attrs.create("header", np.string_(header_text))
for iorder in range(1, 10):
    ecc_data[:, 2*iorder-2] = data_en[:, 4*iorder-2]
    ecc_data[:, 2*iorder-1] = data_en[:, 4*iorder-1]
h5data = hf.create_dataset("ecc_data_en", data=ecc_data.astype('float32'),
                           compression="gzip", compression_opts=9)
h5data.attrs.create("header", np.string_(header_text))

# close the hdf5 file
hf.close()

print("finished.")
