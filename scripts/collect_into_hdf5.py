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
coll_data = np.zeros([nev, 5])
coll_data[:, 0] = data_sn[:, 48]  # b [fm]
coll_data[:, 1] = data_sn[:, 45]  # Npart
coll_data[:, 2] = data_sn[:, 46]  # Ncoll
coll_data[:, 3] = data_sn[:, 47]  # dS/dy
coll_data[:, 4] = data_en[:, 47]  # dE/dy

# save collision data
header_text = "# b[fm]  Npart  Ncoll  tau0*dS/deta_s  tau0*dE/deta_s [GeV]"
h5data = hf.create_dataset("collision_data", data=coll_data.astype('float32'),
                           compression="gzip", compression_opts=9)
h5data.attrs.create("header", np.string_(header_text))

# save <r^n> data
rn_data = np.zeros([nev, 9])
for iorder in range(1, 10):
    rn_data[:, iorder-1] = data_sn[:, 5*iorder-1]
header_text = "# <r^n> [fm^n] (n = 1-9)"
h5data = hf.create_dataset("rn_data_sn", data=rn_data.astype('float32'),
                           compression="gzip", compression_opts=9)
h5data.attrs.create("header", np.string_(header_text))
for iorder in range(1, 10):
    rn_data[:, iorder-1] = data_en[:, 5*iorder-1]
h5data = hf.create_dataset("rn_data_en", data=rn_data.astype('float32'),
                           compression="gzip", compression_opts=9)
h5data.attrs.create("header", np.string_(header_text))

# save eccentricity data
ecc_data = np.zeros([nev, 18])
for iorder in range(1, 10):
    ecc_data[:, 2*iorder-2] = data_sn[:, 5*iorder-3]
    ecc_data[:, 2*iorder-1] = data_sn[:, 5*iorder-2]
header_text = "# ecc_n_real  ecc_n_imag (n = 1-9)"
h5data = hf.create_dataset("ecc_data_sn", data=ecc_data.astype('float32'),
                           compression="gzip", compression_opts=9)
h5data.attrs.create("header", np.string_(header_text))
for iorder in range(1, 10):
    ecc_data[:, 2*iorder-2] = data_en[:, 5*iorder-3]
    ecc_data[:, 2*iorder-1] = data_en[:, 5*iorder-2]
h5data = hf.create_dataset("ecc_data_en", data=ecc_data.astype('float32'),
                           compression="gzip", compression_opts=9)
h5data.attrs.create("header", np.string_(header_text))

# close the hdf5 file
hf.close()

print("finished.")
