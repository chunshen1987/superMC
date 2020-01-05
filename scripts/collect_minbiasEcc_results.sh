#!/usr/bin/env bash

workpath=$1
workpath=${workpath%"/"}
resFolder=`echo ${workpath} | rev | cut -f 1 -d "/" | rev`

echo "collecting results from ${workpath} ..."
mkdir -p ${resFolder}
for ijob in `ls --color=none ${workpath}`
do
    cat ${workpath}/${ijob}/data/en_ecc_eccp_10.dat >> ${resFolder}/en_ecc_eccp_10.dat
    cat ${workpath}/${ijob}/data/sn_ecc_eccp_10.dat >> ${resFolder}/sn_ecc_eccp_10.dat
    cp ${workpath}/${ijob}/test.log ${resFolder}/run_${ijob}.log
    cp ${workpath}/${ijob}/test.err ${resFolder}/run_${ijob}.err
done

./collect_into_hdf5.py ${resFolder}
./centrality_cut_h5.py ${resFolder}/${resFolder}.h5
