#!/bin/bash

if [ "$PS1" ]; then echo -e "This script cannot be sourced. Use \"${BASH_SOURCE[0]}\" instead." ; return ; fi

CONFIGFILE=$1
if [ -z ${CONFIGFILE} ] || [ ! -f ${CONFIGFILE} ]; then
    echo "Fatal! No configuration file found.";
    echo "Use this script as: ${BASH_SOURCE[0]} CONFIGFILE";
    exit 1;
fi
source ${CONFIGFILE}
source PATHS
source EXTERNAL
source ${UTILSDIR}/read_simparams

echo "Number of datasets: ${#SIMPARAMS[@]}"

for PARAMSTR in ${SIMPARAMS[@]}; do
    echo ${PARAMSTR}
    ENDSIM=$(( STARTSIM + NSIM ))

    ALLMETHODS=""
    if [ "${bMatrixEqtl}" = "true" ];  then ALLMETHODS="matrixeqtl ${ALLMETHODS}"; fi
    #if [ "${bTejaas}" = "true" ];      then ALLMETHODS="rr ${ALLMETHODS}"; fi
    #if [ "${bTejaasJPA}" = "true" ];   then ALLMETHODS="jpa ${ALLMETHODS}"; fi

    for THISMETHOD in ${ALLMETHODS}; do
        ${PYTHON36} ${ROCSUMMPY} --startsim ${STARTSIM} \
                                 --endsim ${ENDSIM} \
                                 --method ${THISMETHOD} \
                                 --which ${WHICHPLOTS} \
                                 --srcdir ${OUTDIRUP}/${PARAMSTR} \
                                 --outdir ${OUTDIRUP}/${PARAMSTR}/rocdata
    done

done
