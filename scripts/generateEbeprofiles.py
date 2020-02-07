#! /usr/bin/env python

import sys, shutil
from numpy import *
from os import path, makedirs
import subprocess
import re
from glob import glob

class color:
    """
    define colors in the terminal
    """
    purple = '\033[95m'
    cyan = '\033[96m'
    darkcyan = '\033[36m'
    blue = '\033[94m'
    green = '\033[92m'
    yellow = '\033[93m'
    red = '\033[91m'
    bold = '\033[1m'
    underline = '\033[4m'
    end = '\033[0m'


#dictionary for parameter list in superMC
superMCParameters = {
    'which_mc_model'                :   5,
    'sub_model'                     :   1,
    'Npmin'                         :   2,
    'Npmax'                         :   1000,
    'bmin'                          :   0,
    'bmax'                          :   20,
    'cutdSdy'                       :   1,
    'cutdSdy_lowerBound'            :   551.864,
    'cutdSdy_upperBound'            :   1000000.0,
    'Aproj'                         :   197,
    'Atarg'                         :   197,
    'ecm'                           :   200,
    'finalFactor'                   :   1.0,
    'use_ed'                        :   0,
    'use_sd'                        :   1,
    'alpha'                         :   0.14,
    'lambda'                        :   0.138,
    'operation'                     :   2,
    'cc_fluctuation_model'          :   6,
    'output_TATB'                   :   0,
    'output_rho_binary'             :   0,
    'output_TA'                     :   1,
    'output_rhob'                   :   0,
    'output_spectator_density'      :   1,
    'generate_reaction_plane_avg_profile'  :  0,
    'nev'                           :   1000,
    'average_from_order'            :   2,
    'average_to_order'              :   2,
}

nucleus_name_dict = {
    208: 'Pb',
    197: 'Au',
    238: 'U',
    63: 'Cu',
    1: 'p',
    2: 'd',
    3: 'He3',
}

nucleus_number_dict = {
    'Pb': 208,
    'Au': 197,
    'U': 238,
    'Cu': 63,
    'p': 1,
    'd': 2,
    'He3': 3,
}


def form_assignment_string_from_dict(adict):
    """
        Generate a parameter-equals-value string from the given dictionary. The
        generated string has a leading blank.
    """
    result = ""
    for aparameter in adict.keys():
        result += " {}={}".format(aparameter, adict[aparameter])
    return result


def translate_centrality_cut(centrality_bound, cut_type='total_entropy'):
    """
    translate the centrality boundaries to Npart, dS/dy, b values and update
    the parameter lists for simulations
    """

    centrality_lower_bound = centrality_bound[0]
    centrality_upper_bound = centrality_bound[1]

    if superMCParameters['which_mc_model'] == 5:
        model_name = 'MCGlb'
    elif superMCParameters['which_mc_model'] == 1:
        model_name = 'MCKLN'

    if superMCParameters['cc_fluctuation_model'] != 0:
        multiplicity_fluctuation = 'withMultFluct'
    else:
        multiplicity_fluctuation = 'noMultFluct'

    collision_energy = '%g' % superMCParameters['ecm']

    Aproj = superMCParameters['Aproj']
    Atrag = superMCParameters['Atarg']

    if Aproj == Atrag:  #symmetric collision
        nucleus_name = nucleus_name_dict[Aproj]+nucleus_name_dict[Atrag]
    else:  # asymmetric collision
        nucleus_name = (nucleus_name_dict[min(Aproj, Atrag)]
                        + nucleus_name_dict[max(Aproj, Atrag)])

    centrality_cut_file_name = (
        'iebe_centralityCut_%s_%s_sigmaNN_gauss_d0.9_%s.dat'
        % (cut_type, model_name + nucleus_name + collision_energy,
           multiplicity_fluctuation)
    )

    try:
        centrality_cut_file = loadtxt(
            path.join(path.abspath('./centrality_cut_tables'),
                      centrality_cut_file_name))
    except IOError:
        print("Can not find the centrality cut table for the collision system")
        print(centrality_cut_file_name)
        exit(1)

    lower_idx = (
        centrality_cut_file[:, 0].searchsorted(centrality_lower_bound+1e-30))
    upper_idx = (
        centrality_cut_file[:, 0].searchsorted(centrality_upper_bound))

    cut_value_upper = (
        (centrality_cut_file[lower_idx-1, 1]
         - centrality_cut_file[lower_idx, 1])
        /(centrality_cut_file[lower_idx-1, 0]
          - centrality_cut_file[lower_idx, 0])
        *(centrality_lower_bound - centrality_cut_file[lower_idx-1, 0])
        + centrality_cut_file[lower_idx-1, 1]
    )
    cut_value_low = (
        (centrality_cut_file[upper_idx-1, 1]
         - centrality_cut_file[upper_idx, 1])
        /(centrality_cut_file[upper_idx-1, 0]
          - centrality_cut_file[upper_idx, 0])
        *(centrality_upper_bound - centrality_cut_file[upper_idx-1, 0])
        + centrality_cut_file[upper_idx-1, 1]
    )
    if cut_type == 'total_entropy':
        superMCParameters['cutdSdy'] = 1
        npart_min = min(centrality_cut_file[lower_idx-1:upper_idx+1, 2])
        npart_max = max(centrality_cut_file[lower_idx-1:upper_idx+1, 3])
        b_min = min(centrality_cut_file[lower_idx-1:upper_idx+1, 4])
        b_max = max(centrality_cut_file[lower_idx-1:upper_idx+1, 5])
        superMCParameters['cutdSdy_lowerBound'] = cut_value_low
        superMCParameters['cutdSdy_upperBound'] = cut_value_upper
    elif cut_type == 'Npart':
        superMCParameters['cutdSdy'] = 0
        b_min = min(centrality_cut_file[lower_idx-1:upper_idx+1, 2])
        b_max = max(centrality_cut_file[lower_idx-1:upper_idx+1, 3])
        npart_min = cut_value_low
        npart_max = cut_value_upper
    superMCParameters['Npmax'] = npart_max
    superMCParameters['Npmin'] = npart_min
    superMCParameters['bmax'] = b_max
    superMCParameters['bmin'] = b_min

    #print out information
    print('-'*80)
    print('%s collisions at sqrt{s} = %s A GeV with %s initial conditions'
          % (nucleus_name , collision_energy, model_name))
    print("Centrality : %g - %g"
          % (centrality_lower_bound, centrality_upper_bound) + r"%")
    print('centrality cut on ', cut_type)
    if cut_type == 'total_entropy':
        print('dS/dy :', cut_value_low, '-', cut_value_upper)
    print("Npart: ", npart_min, '-', npart_max)
    print("b: ", b_min, '-', b_max, ' fm')
    print('-'*80)
    return


def update_superMC_dict(model, ecm, collsys, nev):
    """
    update the superMCParameters dictionary with users input settings
    """
    superMCParameters['nev'] = nev
    if model == 'MCGlb':
        superMCParameters['which_mc_model'] = 5
        superMCParameters['sub_model'] = 1
    elif model == 'MCKLN':
        superMCParameters['which_mc_model'] = 1
        superMCParameters['sub_model'] = 7
        superMCParameters['cc_fluctuation_model'] = 0
    else:
        print(sys.argv[0], ': invalid initial model type', model)
        print_help_message()
        sys.exit(1)

    superMCParameters['ecm'] = ecm
    if ecm == 2760:
        if model == 'MCGlb':
            superMCParameters['alpha'] = 0.118
        elif model == 'MCKLN':
            superMCParameters['lambda'] = 0.138
    if ecm <= 200:
        if model == 'MCGlb':
            superMCParameters['alpha'] = 0.14
        elif model == 'MCKLN':
            superMCParameters['lambda'] = 0.218
    superMCParameters['Aproj'] = nucleus_number_dict[collsys[0]]
    superMCParameters['Atarg'] = nucleus_number_dict[collsys[1]]

    # for checking
    #for x in superMCParameters.keys():
    #    print x + ': ' + str(superMCParameters[x])

    return

def generateEbeprofiles(output_path, centrality_bounds, 
                        cut_type='total_entropy'):
    runRecord = open('./runRecord.dat', 'a')
    errRecord = open('./errRecord.dat', 'a')
    if not path.exists(output_path):
        makedirs(output_path)
    translate_centrality_cut(centrality_bounds, cut_type)
    cen_string = '%g-%g' %(centrality_bounds[0], centrality_bounds[1])
    option = form_assignment_string_from_dict(superMCParameters)
    cmd = './superMC.e' + option
    superMC_folder = path.abspath('./')
    print(cmd)
    runRecord.write(cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=runRecord,
                         stderr=errRecord, cwd=superMC_folder)
    p.wait()

    # save files
    store_folder = output_path
    filelist = glob(path.join(superMC_folder, "data", "*"))
    for ifile, filepath in enumerate(filelist):
        filename = filepath.split("/")[-1]
        if "block" in filename:
            newfilename = re.sub("event", "C{}_event".format(cen_string),
                                 filename)
            newfilename = re.sub("_block", "".format(cen_string), newfilename)
            shutil.move(filepath, path.join(store_folder, newfilename))
    shutil.move('./runRecord.dat', path.join(store_folder, 'runRecord.dat'))
    shutil.move('./errRecord.dat', path.join(store_folder, 'errRecord.dat'))

def print_help_message():
    print("Usage : ")
    print(color.bold
          + "{} -ecm ecm ".format(sys.argv[0])
          + "-cen cen_bounds"
          + "[-model model -collision_system collsys -cut_type cut_type]"
          + color.end)
    print("Usage of {} command line arguments: ".format(sys.argv[0]))
    print(color.bold + "-cen" + color.end
          + "   centrality bounds(%): "
          + color.purple + "20-30" + color.end)
    print(color.bold + "-ecm" + color.end
          + "   collision energy (GeV): "
          + color.purple + "7.7, 11.5, 19.6, 27, 39, 62.4, 200, 2760, 5500"
          + color.end)
    print(color.bold + "-cut_type" + color.end
          + "   centrality cut type: "
          + color.purple + color.bold + "total_entropy[default]" + color.end
          + color.purple + ", Npart" + color.end)
    print(color.bold + "-model" + color.end + " initial condition model: "
          + color.purple + color.bold + " MCGlb[default]" + color.end
          + color.purple + ", MCKLN" + color.end)
    print(color.bold + "-collision_system" + color.end
          + " type of collision system: "
          + color.purple + color.bold + " Pb+Pb[default]" + color.end
          + color.purple + ", Au+Au, Cu+Au, U+U, p+Pb, p+Au, d+Au, He3+Au"
          + color.end)


if __name__ == "__main__":
    # set default values
    model = 'MCGlb'
    cut_type = 'total_entropy'
    #cut_type = 'Npart'
    collsys = 'Au+Au'.split('+')
    output_path = path.abspath('./RESULTS/')
    nev = 100

    while len(sys.argv) > 1:
        option = sys.argv[1]
        del sys.argv[1]
        if option == '-model':
            model = str(sys.argv[1])
            del sys.argv[1]
        elif option == '-collision_system':
            collsys = str(sys.argv[1]).split('+')
            del sys.argv[1]
        elif option == "-nev":
            nev = int(sys.argv[1])
            del sys.argv[1]
        elif option == '-cut_type':
            cut_type = str(sys.argv[1])
            del sys.argv[1]
            if cut_type not in ['total_entropy', 'Npart']:
                print(sys.argv[0], ': invalid centrality cut type', cut_type)
                print_help_message()
                sys.exit(1)
        elif option == '-cen':
            centrality_bounds = [float(istr) for istr in str(sys.argv[1]).split('-')]
            del sys.argv[1]
        elif option == '-ecm':
            ecm = float(sys.argv[1])
            del sys.argv[1]
        elif option == '-output':
            folder = float(sys.argv[1])
            output_path = path.join(path.abspath('./'), folder)
            del sys.argv[1]
        elif option == '-h':
            print_help_message()
            sys.exit(0)
        else:
            print(sys.argv[0], ': invalid option ', option)
            print_help_message()
            sys.exit(1)

    try:
        update_superMC_dict(model, ecm, collsys, nev)
        generateEbeprofiles(output_path, centrality_bounds, cut_type)
    except NameError:
        print_help_message()
        sys.exit(1)
