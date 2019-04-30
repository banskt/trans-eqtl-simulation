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
        for NULL in ${TEJAAS_NULL}; do
            if [ ${NULL} == "perm" ]; then TEJAAS_SIGMA_BETA=${TEJAAS_SIGMA_BETA_PERM}; fi
            if [ ${NULL} == "maf" ]; then TEJAAS_SIGMA_BETA=${TEJAAS_SIGMA_BETA_MAF}; fi
            for SBETA in ${TEJAAS_SIGMA_BETA}; do
                for NPCA in ${TEJAAS_PRINCIPAL_COMPONENTS}; do
                    for NPEER in ${NPEERCORR}; do
                        if [ ! "${NPEER}" = "None" ]; then 
                            if [ "${bTejaas}" = "true" ];    then tejaas_chunk_reduce "${OUTDIR_SIM}/tejaas/${NULL}null_sb${SBETA}_${NPCA}pc/npeer${NPEER}" $NCHUNK; fi
                            if [ "${bTjsRandom}" = "true" ]; then tejaas_chunk_reduce "${OUTDIR_SIM}/tejaas_rand/${NULL}null_sb${SBETA}_${NPCA}pc/npeer${NPEER}" $NCHUNK; fi
                        else
                            if [ "${bTejaas}" = "true" ];    then tejaas_chunk_reduce "${OUTDIR_SIM}/tejaas/${NULL}null_sb${SBETA}_${NPCA}pc" $NCHUNK; fi
                            if [ "${bTjsRandom}" = "true" ]; then tejaas_chunk_reduce "${OUTDIR_SIM}/tejaas_rand/${NULL}null_sb${SBETA}_${NPCA}pc" $NCHUNK; fi
                        fi
                    done
                done
            done
        done
    done

done
