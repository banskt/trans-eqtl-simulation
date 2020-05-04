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

CONFIGSTRING=$( basename ${CONFIGFILE} )
CONFIGSTRING=${CONFIGSTRING#"CONFIG."}
echo "Number of datasets: ${#SIMPARAMS[@]}"

RANDSTRING=$( cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 4 | head -n 1 )

THISJOBNAME="rocdata_${CONFIGSTRING}_${RANDSTRING}"
SPECIFIC_JOBSUBDIR=${JOBSUBDIR}/rocdata
THISJOBFILE="${SPECIFIC_JOBSUBDIR}/${THISJOBNAME}.sbatch"
if [ ! -d ${SPECIFIC_JOBSUBDIR} ]; then mkdir -p ${SPECIFIC_JOBSUBDIR}; fi
sed -e "s|_JOB_NAME|${THISJOBNAME}|g;" ${MASTER_BSUBDIR}/rocdata.bsub > ${THISJOBFILE}
echo " " >> ${THISJOBFILE}

for PARAMSTR in ${SIMPARAMS[@]}; do

    echo ${PARAMSTR}


    ENDSIM=$(( STARTSIM + NSIM ))
    NSNP=$( echo ${PARAMSTR} | cut -d'_' -f1 )
    NTOP1=$(( NSNP / 10 ))
    NTOP2=$(( NSNP / 20 ))
    NTOP3=$(( NSNP / 100 ))

    SCRIPTBEGIN="${PYTHON36} ${ROCSUMMPY} --startsim ${STARTSIM} --endsim ${ENDSIM} --which ${WHICHPLOTS} --ntop ${NTOP1} ${NTOP2} ${NTOP3} --srcdir ${OUTDIRUP}/${PARAMSTR} --outdir ${OUTDIRUP}/${PARAMSTR}/rocdata"

    for NPEER in ${MEQTL_NPEER}; do
        for PRCC in ${MATRIXEQTL_PREPROC}; do
            if [ "${bMatrixEqtl}" = "true" ]; then
                echo "${SCRIPTBEGIN} --method matrixeqtl --mdir matrixeqtl/${PRCC}/peer${NPEER} --outprefix matrixeqtl_${PRCC}_peer${NPEER}" >> ${THISJOBFILE}
            fi
            if [ "${bMEqtlRandom}" = "true" ]; then
                echo "${SCRIPTBEGIN} --method matrixeqtl --mdir matrixeqtl_rand/${PRCC}/peer${NPEER} --outprefix matrixeqtl_rand_${PRCC}_peer${NPEER}" >> ${THISJOBFILE}
            fi
        done
    done

    for NPEER in ${JPA_NPEER}; do
        for PRIDX in ${!JPA_PREPROC[*]}; do
            PRCC=${JPA_PREPROC[PRIDX]}
            KNN_NBR=${JPA_KNN[PRIDX]}
            if [ ! "${KNN_NBR}" = "0" ]; then PRCC="${PRCC}_knn${KNN_NBR}"; fi
            if [ "${bJPA}" = "true" ]; then
                echo "${SCRIPTBEGIN} --method jpa --mdir jpa/${PRCC}/peer${NPEER} --outprefix jpa_${PRCC}_peer${NPEER}" >> ${THISJOBFILE}
            fi
            if [ "${bJPARandom}" = "true" ]; then
                echo "${SCRIPTBEGIN} --method jpa --mdir jpa_rand/${PRCC}/peer${NPEER} --outprefix jpa_rand_${PRCC}_peer${NPEER}" >> ${THISJOBFILE}
            fi
        done
    done

    for NPEER in ${TEJAAS_NPEER}; do
        for PRIDX in ${!TEJAAS_PREPROC[*]}; do
            PRCC=${TEJAAS_PREPROC[PRIDX]}
            KNN_NBR=${TEJAAS_KNN[PRIDX]}
            for NULL in ${TEJAAS_NULL}; do
                if [ ${NULL} == "perm" ]; then __SIGMA_BETA=${TEJAAS_SIGMA_BETA_PERM}; fi
                if [ ${NULL} == "maf" ];  then __SIGMA_BETA=${TEJAAS_SIGMA_BETA_MAF}; fi
                for SBETA in ${__SIGMA_BETA}; do
                    METHOD_VARIANT="${NULL}null_sb${SBETA}"
                    if [ ! "${TEJAAS_CISMASK}" = "true" ]; then METHOD_VARIANT="${METHOD_VARIANT}_nomask"; fi
                    PRCCMOD=${PRCC}
                    if [ ! "${KNN_NBR}" = 0 ]; then PRCCMOD="${PRCCMOD}_knn${KNN_NBR}"; fi
                    if [ "${bTejaas}" = "true" ]; then
                        echo "${SCRIPTBEGIN} --method rr --mdir tejaas/${METHOD_VARIANT}/${PRCCMOD}/peer${NPEER} --outprefix tejaas_${METHOD_VARIANT}_${PRCCMOD}_peer${NPEER}" >> ${THISJOBFILE}
                    fi
                    if [ "${bTejaasRandom}" = "true" ]; then
                        echo "${SCRIPTBEGIN} --method rr --mdir tejaas_rand/${METHOD_VARIANT}/${PRCCMOD}/peer${NPEER} --outprefix tejaas_rand_${METHOD_VARIANT}_${PRCCMOD}_peer${NPEER}" >> ${THISJOBFILE}
                    fi
                done
            done
        done
    done
done
submit_job ${SPECIFIC_JOBSUBDIR} ${THISJOBNAME} None
