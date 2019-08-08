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
source ${UTILSDIR}/submit_job

echo "Number of datasets: ${#SIMPARAMS[@]}"

RANDSTRING=$( cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1 )


for PARAMSTR in ${SIMPARAMS[@]}; do

    echo ${PARAMSTR}

    THISJOBNAME="rocdata_${RANDSTRING}_${PARAMSTR}"
    SPECIFIC_JOBSUBDIR=${JOBSUBDIR}/${PARAMSTR}/rocdata
    THISJOBFILE="${SPECIFIC_JOBSUBDIR}/${THISJOBNAME}.sh"
    if [ -d ${SPECIFIC_JOBSUBDIR} ]; then rm -rf ${SPECIFIC_JOBSUBDIR}; fi; mkdir -p ${SPECIFIC_JOBSUBDIR}
    sed -e "s|_JOB_NAME|${THISJOBNAME}|g;" ${MASTER_BSUBDIR}/rocdata.bsub > ${THISJOBFILE}
    echo " " >> ${THISJOBFILE}
    echo " " >> ${THISJOBFILE}

    ENDSIM=$(( STARTSIM + NSIM ))
    NSNP=$( echo ${PARAMSTR} | cut -d'_' -f1 )
    NTOP1=$(( NSNP / 10 ))
    NTOP2=$(( NSNP / 20 ))
    NTOP3=$(( NSNP / 100 ))

    SCRIPTBEGIN="${PYTHON36} ${ROCSUMMPY} --startsim ${STARTSIM} --endsim ${ENDSIM} --which ${WHICHPLOTS} --ntop ${NTOP1} ${NTOP2} ${NTOP3} --srcdir ${OUTDIRUP}/${PARAMSTR} --outdir ${OUTDIRUP}/${PARAMSTR}/rocdata"

    for NPEER in ${NPEERCORR}; do
        for PRCC in ${MATRIXEQTL_PREPROC}; do
            if [ "${bMatrixEqtl}" = "true" ]; then
                echo "${SCRIPTBEGIN} --method matrixeqtl --mdir matrixeqtl/${PRCC}/peer${NPEER} --outprefix matrixeqtl_${PRCC}_peer${NPEER}" >> ${THISJOBFILE}
            fi
            if [ "${bMEqtlRandom}" = "true" ]; then
                echo "${SCRIPTBEGIN} --method matrixeqtl --mdir matrixeqtl_rand/${PRCC}/peer${NPEER} --outprefix matrixeqtl_rand_${PRCC}_peer${NPEER}" >> ${THISJOBFILE}
            fi
        done
    done

    for NPEER in ${TEJAAS_NPEER}; do
        for PRIDX in ${!TEJAAS_PREPROC[*]}; do
            PRCC=${TEJAAS_PREPROC[PRIDX]}
            KNN=${TEJAAS_KNN[PRIDX]}
            for NULL in ${TEJAAS_NULL}; do
                if [ ${NULL} == "perm" ]; then __SIGMA_BETA=${TEJAAS_SIGMA_BETA_PERM}; fi
                if [ ${NULL} == "maf" ];  then __SIGMA_BETA=${TEJAAS_SIGMA_BETA_MAF}; fi
                for SBETA in ${__SIGMA_BETA}; do
                    METHOD_VARIANT="${NULL}null_sb${SBETA}"
                    if [ ! "${TEJAAS_CISMASK}" = "true" ]; then METHOD_VARIANT="${METHOD_VARIANT}_nomask"; fi
                    if [ "${KNN}" = "true" ]; then METHOD_VARIANT="${METHOD_VARIANT}_knn"; fi
                    if [ "${bTejaas}" = "true" ]; then
                        echo "${SCRIPTBEGIN} --method rr --mdir tejaas/${METHOD_VARIANT}/${PRCC}/peer${NPEER} --outprefix tejaas_${METHOD_VARIANT}_${PRCC}_peer${NPEER}" >> ${THISJOBFILE}
                    fi
                    if [ "${bTejaasRandom}" = "true" ]; then
                        echo "${SCRIPTBEGIN} --method rr --mdir tejaas_rand/${METHOD_VARIANT}/${PRCC}/peer${NPEER} --outprefix tejaas_rand_${METHOD_VARIANT}_${PRCC}_peer${NPEER}" >> ${THISJOBFILE}
                    fi
                done
            done
        done
    done

    #submit_job ${SPECIFIC_JOBSUBDIR} ${THISJOBNAME} None

done
