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
    #echo ${PARAMSTR}
    ENDSIM=$(( STARTSIM + NSIM ))
    NSNP=$( echo ${PARAMSTR} | cut -d'_' -f1 )
    NTOP1=$(( NSNP / 10 ))
    NTOP2=$(( NSNP / 20 ))
    NTOP3=$(( NSNP / 100 ))

    ALLMETHODS=""
    #if [ "${bMatrixEqtl}" = "true" ];  then ALLMETHODS="matrixeqtl ${ALLMETHODS}"; fi
    if [ "${bTejaas}" = "true" ];      then ALLMETHODS="rr ${ALLMETHODS}"; fi
    #if [ "${bTjsRandom}" = "true" ];   then ALLMETHODS="rr_rand ${ALLMETHODS}"; fi
    #if [ "${bTejaasJPA}" = "true" ];   then ALLMETHODS="jpa ${ALLMETHODS}"; fi

    for THISMETHOD in ${ALLMETHODS}; do
        for NPEER in ${NPEERCORR}; do
            if [ ! "${NPEER}" = "None" ]; then PEERFLAG="--npeer ${NPEER}"; fi
            SIGBETALOOP="None"; if [ "${THISMETHOD}" = "rr" ]; then SIGBETALOOP="${TEJAAS_SIGMA_BETA_PERM}"; fi
            #SIGBETALOOP="None"; if [ "${THISMETHOD}" = "rr_rand" ]; then SIGBETALOOP="${TEJAAS_SIGMA_BETA_PERM}"; fi
            for SBETA in ${SIGBETALOOP}; do
                echo "${PYTHON36} ${ROCSUMMPY} --startsim ${STARTSIM} --endsim ${ENDSIM} --method ${THISMETHOD} --which ${WHICHPLOTS} --sbeta ${SBETA} --ntop ${NTOP1} ${NTOP2} ${NTOP3} --srcdir ${OUTDIRUP}/${PARAMSTR} --outdir ${OUTDIRUP}/${PARAMSTR}/rocdata ${PEERFLAG}"
                ${PYTHON36} ${ROCSUMMPY} --startsim ${STARTSIM} \
                                         --endsim ${ENDSIM} \
                                         --method ${THISMETHOD} \
                                         --which ${WHICHPLOTS} \
                                         --sbeta ${SBETA} \
                                         --ntop ${NTOP1} ${NTOP2} ${NTOP3} \
                                         --srcdir ${OUTDIRUP}/${PARAMSTR} \
                                         --outdir ${OUTDIRUP}/${PARAMSTR}/rocdata ${PEERFLAG}
                #echo "$PARAMSTR	$THISMETHOD	$SBETA"
            done
        done
    done

done
