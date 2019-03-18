#!/bin/bash

CURDIR=$( pwd )
while read ERRFILE; do
    #echo $ERRFILE
    JOBSUBDIR=$( dirname $ERRFILE )
    JOBNAME=$( basename $ERRFILE .err_8core )
    echo $JOBSUBDIR
    cd $JOBSUBDIR
        #mv ${JOBNAME}.err ${JOBNAME}.err_8core
        #mv ${JOBNAME}.out ${JOBNAME}.out_8core
        #rm -rf ${JOBNAME}.err ${JOBNAME}.out
        sed -e "s|-n\ 8|-n\ 16|g" ${JOBNAME}.bsub > ${JOBNAME}_16core.bsub
        bsub < ${JOBNAME}_16core.bsub
    cd $CURDIR
done < jpa_error_files.txt
