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
source ${UTILSDIR}/tejaas_chunk_reduce.new

echo "Number of datasets: ${#SIMPARAMS[@]}"

for PARAMSTR in ${SIMPARAMS[@]}; do
    ENDSIM=$(( STARTSIM + NSIM ))

    for (( SIM=$STARTSIM; SIM<$ENDSIM; SIM++ )); do
        SIMINDEX=`echo $SIM | awk '{printf "%03d", $1}'`
        OUTDIR_SIM="${OUTDIRUP}/${PARAMSTR}/sim${SIMINDEX}"
        SIMGTFILE="${OUTDIR_SIM}/input/genotype.vcf.gz"
        NCHUNK=$( expected_nchunk ${SIMGTFILE} ${MAX_NSNP_PERJOB} )
        echo "Simulation ${SIMINDEX} --> ${NCHUNK} nchunks"

        for PRIDX in ${!TEJAAS_PREPROC[*]}; do

            PRCC=${TEJAAS_PREPROC[PRIDX]}
            KNN=${TEJAAS_KNN[PRIDX]}

            for NULL in ${TEJAAS_NULL}; do
                if [ ${NULL} == "perm" ]; then __SIGMA_BETA=${TEJAAS_SIGMA_BETA_PERM}; fi
                if [ ${NULL} == "maf" ]; then  __SIGMA_BETA=${TEJAAS_SIGMA_BETA_MAF}; fi
                for SBETA in ${__SIGMA_BETA}; do
                    METHOD_VARIANT="${NULL}null_sb${SBETA}"
                    if [ "${KNN}" = "true" ]; then
                        METHOD_VARIANT="${METHOD_VARIANT}_knn"
                    fi
                    for NPEER in ${TEJAAS_NPEER}; do
                            if [ "${bTejaas}" = "true" ];    then tejaas_chunk_reduce "${OUTDIR_SIM}/tejaas/${METHOD_VARIANT}/${PRCC}/peer${NPEER}" $NCHUNK; fi
                            if [ "${bTjsRandom}" = "true" ]; then tejaas_chunk_reduce "${OUTDIR_SIM}/tejaas_rand/${METHOD_VARIANT}/${PRCC}/peer${NPEER}" $NCHUNK; fi
                    done
                done
            done
        done
    done
done
