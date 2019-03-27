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
source ${UTILSDIR}/read_simparams
source ${UTILSDIR}/tejaas_chunk_reduce

echo "Number of datasets: ${#SIMPARAMS[@]}"

for PARAMSTR in ${SIMPARAMS[@]}; do
    ENDSIM=$(( STARTSIM + NSIM ))

    for (( SIM=$STARTSIM; SIM<$ENDSIM; SIM++ )); do
        SIMINDEX=`echo $SIM | awk '{printf "%03d", $1}'`
        OUTDIR_SIM="${OUTDIRUP}/${PARAMSTR}/sim${SIMINDEX}"
        echo "Simulation ${SIMINDEX}"
        ## if [ "${bTejaasJPA}" = "true" ]; then tejaas_chunk_reduce "${OUTDIR_SIM}/tejaas/jpa"; fi
        ## if [ "${bJPARandom}" = "true" ]; then tejaas_chunk_reduce "${OUTDIR_SIM}/tejaas_rand/jpa"; fi
        for NULL in ${TEJAAS_NULL}; do
            if [ ${NULL} == "perm" ]; then TEJAAS_SIGMA_BETA=${TEJAAS_SIGMA_BETA_PERM}; fi
            if [ ${NULL} == "maf" ]; then TEJAAS_SIGMA_BETA=${TEJAAS_SIGMA_BETA_MAF}; fi
            for SBETA in ${TEJAAS_SIGMA_BETA}; do
                for NPEER in ${NPEERCORR}; do
                    if [ ! "${NPEER}" = "None" ]; then 
                        if [ "${bTejaas}" = "true" ];    then tejaas_chunk_reduce "${OUTDIR_SIM}/tejaas/${NULL}null_sb${SBETA}/npeer${NPEER}"; fi
                        if [ "${bTjsRandom}" = "true" ]; then tejaas_chunk_reduce "${OUTDIR_SIM}/tejaas_rand/${NULL}null_sb${SBETA}/npeer${NPEER}"; fi
                    else
                        if [ "${bTejaas}" = "true" ];    then tejaas_chunk_reduce "${OUTDIR_SIM}/tejaas/${NULL}null_sb${SBETA}"; fi
                        if [ "${bTjsRandom}" = "true" ]; then tejaas_chunk_reduce "${OUTDIR_SIM}/tejaas_rand/${NULL}null_sb${SBETA}"; fi
                    fi
                done
            done
        done
    done

done
