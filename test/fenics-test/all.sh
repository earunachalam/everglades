#!/bin/bash --login

TASKNAME=fish-pigmentation
BASEDIR=$(pwd)
DATADIR=dat/${TASKNAME}/
IMGDIR=img/${TASKNAME}/

# delete previously created csv files
(cd ${DATADIR} && rm -f *)

# run in fenics
python 2d.py

# export pvd files to csv using paraview
(cd ${DATADIR} && ${BASEDIR}/pvd2csv.py)

# rename files to eliminate periods (generated automatically by paraview export)
(cd ${DATADIR} && prename 's/(u_[0-9]+)\.([0-9]+\.csv)/$1_$2/' *)

# create images
(cd ${IMGDIR} && rm *.pdf)
NTSTEP=$(bc <<< "$(ls dat/fish-pigmentation/*.csv | wc -l)/2")
INTTSTEP=5
(cd plt && ./plot_2d_colored.py ${DATADIR} ${IMGDIR} ${NTSTEP} ${INTTSTEP})

echo Completed.
