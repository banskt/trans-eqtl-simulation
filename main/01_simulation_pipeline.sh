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
source ${UTILSDIR}/submit_job
source ${UTILSDIR}/add_deps
source ${UTILSDIR}/read_simparams

echo "Number of datasets: ${#SIMPARAMS[@]}"

# Run simulation for each set of parameters
for PARAMSTR in ${SIMPARAMS[@]}; do
    echo "${PARAMSTR}"
    ENDSIM=$(( STARTSIM + NSIM ))
    for (( SIM=$STARTSIM; SIM<$ENDSIM; SIM++ )); do


        SIMINDEX=`echo $SIM | awk '{printf "%03d", $1}'`
        RANDSTRING=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1`
        SHUFFLE=false   # used for controlling shuffling
        SUBMITTED_JOBIDS="" # used for controlling jobid reporting

        ## control job dependencies
        GENDATA_JOBDEPS="None"
        PREPROC_JOBDEPS="None"
        MATRIXEQTL_JOBDEPS="None"
        TEJAAS_JOBDEPS="None"

        JOBSUBDIR_SIM="${JOBSUBDIR}/${PARAMSTR}/sim${SIMINDEX}"
        OUTDIR_SIM="${OUTDIRUP}/${PARAMSTR}/sim${SIMINDEX}"
        SIMGTFILE="${OUTDIR_SIM}/input/genotype.vcf.gz"
        SIMGXFILE="${OUTDIR_SIM}/input/expression.txt"
        SIMCFFILE="${OUTDIR_SIM}/input/expression.cf"
        GXPROCFILE="${OUTDIR_SIM}/input/gx.txt"

        echo "  sim${SIMINDEX}:"

        if [ "${bGenerateData}" = "true" ]; then 
            source ${UTILSDIR}/generate_data
            echo "    ${GENDATA_JOBID} > Generate data."
        fi

        if [ "${bPreprocessData}" = "true" ]; then 
            source ${UTILSDIR}/preprocess_data
            echo "    ${PREPROC_JOBID} > Preprocess data: ${GENDATA_JOBDEPS}"
        fi

        for NPEER in ${NPEERCORR}; do
            if [ "${bMatrixEqtl}" = "true" ];   then source ${UTILSDIR}/matrix_eqtl; fi
            if [ "${bMEqtlRandom}" = "true" ];  then SHUFFLE=true; source ${UTILSDIR}/matrix_eqtl; SHUFFLE=false; fi
            #if [ "${bJPA}" = "true" ];          then source ${UTILSDIR}/jpa; fi
            #if [ "${bJPARandom}" = "true" ];    then SHUFFLE=true; source ${UTILSDIR}/jpa; SHUFFLE=false; fi
        done

        for NPEER in ${TEJAAS_NPEER}; do
            if [ "${bTejaas}" = "true" ];       then source ${UTILSDIR}/tejaas; fi
            if [ "${bTejaasRandom}" = "true" ]; then SHUFFLE=true; source ${UTILSDIR}/tejaas; SHUFFLE=false; fi
        done
    done
done
