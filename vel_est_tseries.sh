#!/bin/bash

RUN_VEL_EST=${HOME}/src/velocity_imaging/run_vel_est


START=0
STOP=15
STID=0

LONG_OPTIONS="stid,start_beam:,stop_beam:"

OPTIONS=$(getopt -l "$LONG_OPTIONS" -- "$@")

while true; do
    case "$1" in
	--stid)
	    STID=1
	    shift
	    ;;
	--start_beam)
	    START="$2"
	    shift 2
	    ;;
	--stop_beam)
	    STOP="$2"
	    shift 2
	    ;;
	--)
	    shift
	    break
	    ;;
	*) break ;;
    esac
done


DATE=$1
RADAR=$2
DIR=$3
RCHAR=$4
START_T=$5
END_T=$6

flist=`ls ${DIR}/${DATE}*.iraw.${RCHAR}`

echo ${flist}

for f in ${flist}
do
    echo $f
    strs=$(echo $f | tr "/" "\n")
    for st in $strs
    do
	fdate=${st:0:12}
	ftime=${st:8:12}
    done
    if  [ "${ftime}" \> "${START_T}" ] && [ "${ftime}" \< "${END_T}" ]; then

	if [ $STID -eq 0 ]
	then
	    ${RUN_VEL_EST} --start_beam ${START} --stop_beam ${STOP} ${f} ${RADAR} ${DIR} ${RCHAR}
	else
	    ${RUN_VEL_EST} --start_beam ${START} --stop_beam ${STOP} --stid ${f} ${RADAR} ${DIR} ${RCHAR}
	fi		
	v_image_merge ${fdate} ${RADAR}
    fi
done	 

